// =============================================================================
// sparsity_scanner.sv   (v3 — edge-list ingest, Gset scale)
//
// ─────────────────────────────────────────────────────────────────────────────
// WHY THIS CHANGED FROM THE PROTOTYPE
// ─────────────────────────────────────────────────────────────────────────────
// The prototype scanned a DENSE N×N grid (N*N reads) looking for non-zero
// J[i][j].  At Gset scale that is N*N = 400,000,000 reads for N=20,000 —
// completely impractical, and pointless since we now store only edges.
//
// This version instead walks the EDGE LIST ONCE (num_edges reads, typically
// a few tens of thousands even for the largest Gset instance) and builds the
// CSR store directly.  Since Gset graphs are undirected, each edge (i,j,w)
// contributes to BOTH row i (entry pointing to j) and row j (entry pointing
// to i) of the symmetric J matrix.
//
// THE HARD PART — building row_ptr without knowing row order in advance
// -----------------------------------------------------------------------
// Edges arrive in arbitrary order (whatever order the Gset file lists them).
// CSR requires each row's entries to be CONTIGUOUS in the entry array, but
// we don't know how many entries row i will have until we've seen every
// edge.  This module solves it in TWO PASSES over the edge list:
//
//   PASS 1 (COUNT)   : walk all edges, increment degree[i] and degree[j]
//                       for each edge (no CSR writes yet)
//   FINALIZE          : convert degree[] into row_ptr[] start offsets via
//                       a running prefix sum (row_ptr[i].start = sum of
//                       degree[0..i-1]); count field = degree[i]
//   PASS 2 (PLACE)    : walk all edges again; for edge (i,j,w) write the
//                       entry for row i at  (row_ptr[i].start + fill[i])
//                       and increment fill[i]; same for row j with (j,i,w).
//                       fill[] is a separate write-cursor array, reset to 0
//                       after FINALIZE, incremented as each row fills up.
//
// This is the standard two-pass CSR construction algorithm (counting sort
// style) and requires only O(N) extra storage (degree/fill arrays) and
// O(2 * num_edges) total reads — entirely practical even at Gset's largest
// scale (82,918 directed entries for G63/G64).
//
// degree[] and fill[] are stored as internal block RAM (N entries each,
// ROW_COUNT_W bits wide) — sized for N up to dsb_pkg::N.
//
// Timing
// ------
//  main_memory has 1-cycle read latency (registered output) — same pipeline
//  shape as the prototype scanner: READ -> WAIT -> (eval/write) -> next.
// =============================================================================

`timescale 1ns/1ps
module sparsity_scanner
    import dsb_pkg::*;
    import dsb_sparse_pkg::*;
#(
    parameter int unsigned N_P          = N,
    parameter int unsigned IC_BITS_P    = IC_BITS,
    parameter int unsigned MAIN_ADDR_WP = MAIN_ADDR_W,
    parameter int unsigned MAIN_WORD_WP = MAIN_WORD_W,
    parameter int unsigned COL_IDX_W_P  = COL_IDX_W,
    parameter int unsigned CSR_ADDR_W_P = CSR_ADDR_W,
    parameter int unsigned CSR_ENTRY_WP = CSR_ENTRY_W,
    parameter int unsigned ROW_PTR_W_P  = ROW_PTR_W,
    parameter int unsigned ROW_COUNT_WP = ROW_COUNT_W
)(
    input  logic clk,
    input  logic rst_n,
    input  logic start,                 // pulse high for one cycle to begin

    // ── Number of edges actually present in main_memory's edge-list region ───
    // (Gset instances vary widely in edge count; this lets one build support
    //  every instance without recompiling MAX_EDGES per file.)
    input  logic [CSR_ADDR_W_P-1:0]  num_edges,

    // ── Active oscillator count for THIS run (<= N_P) ─────────────────────────
    // Allows running an 800-node Gset instance inside an N_P=20000 build
    // without wasting time finalizing unused rows.
    input  logic [COL_IDX_W_P-1:0]   active_n,

    // ── main_memory read port ─────────────────────────────────────────────────
    output logic                     mm_rd_en,
    output logic [MAIN_ADDR_WP-1:0]  mm_rd_addr,
    input  logic [MAIN_WORD_WP-1:0]  mm_rd_data,
    input  logic                     mm_rd_valid,

    // ── csr_index_store write ports ───────────────────────────────────────────
    output logic                     rp_wr_en,
    output logic [COL_IDX_W_P-1:0]   rp_wr_row,
    output logic [ROW_PTR_W_P-1:0]   rp_wr_data,

    output logic                     en_wr_en,
    output logic [CSR_ADDR_W_P-1:0]  en_wr_addr,
    output logic [CSR_ENTRY_WP-1:0]  en_wr_data,

    // ── Status ────────────────────────────────────────────────────────────────
    output logic                     done,
    output logic [CSR_ADDR_W_P-1:0]  nnz_total    // total directed CSR entries
);

    // ── Internal degree / fill / start arrays (one entry per oscillator) ──────
    logic [ROW_COUNT_WP-1:0] degree [0:N_P-1];   // pass 1 output
    logic [CSR_ADDR_W_P-1:0] rstart [0:N_P-1];   // finalize output (prefix sum)
    logic [ROW_COUNT_WP-1:0] fill   [0:N_P-1];   // pass 2 write cursor

    // ── FSM ───────────────────────────────────────────────────────────────────
    typedef enum logic [3:0] {
        ST_IDLE,
        ST_CLEAR,         // zero degree[]/fill[] for active_n rows
        ST_P1_READ, ST_P1_WAIT, ST_P1_ACC,         // pass 1: count degrees
        ST_FIN_INIT,                                // begin prefix sum
        ST_FIN_STEP,                                // one prefix-sum step per row
        ST_FIN_WRROWPTR,                            // write row_ptr[i] to CSR store
        ST_P2_READ, ST_P2_WAIT, ST_P2_PLACE_FWD, ST_P2_PLACE_REV, // pass 2: place entries
        ST_DONE
    } scan2_state_t;

    scan2_state_t state;

    logic [CSR_ADDR_W_P-1:0] edge_idx;     // current edge index (0..num_edges-1)
    logic [COL_IDX_W_P-1:0]  clear_idx;    // row index used during CLEAR / FINALIZE
    logic [CSR_ADDR_W_P-1:0] running_sum;  // running prefix sum during FINALIZE

    // Latched edge fields (unpacked from mm_rd_data)
    logic [COL_IDX_W_P-1:0] e_row, e_col;
    logic [MAIN_ADDR_WP-1:0] e_addr;       // address of this edge word (= edge_idx)

    // ── Edge word unpack ──────────────────────────────────────────────────────
    // word = {row[COL_IDX_W-1:0], col[COL_IDX_W-1:0], weight[IC_BITS-1:0]}
    function automatic logic [COL_IDX_W_P-1:0] unpack_row(logic [MAIN_WORD_WP-1:0] w);
        return w[MAIN_WORD_WP-1 -: COL_IDX_W_P];
    endfunction
    function automatic logic [COL_IDX_W_P-1:0] unpack_col(logic [MAIN_WORD_WP-1:0] w);
        return w[IC_BITS_P +: COL_IDX_W_P];
    endfunction

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state        <= ST_IDLE;
            edge_idx     <= '0;
            clear_idx    <= '0;
            running_sum  <= '0;
            done         <= 1'b0;
            nnz_total    <= '0;
            mm_rd_en     <= 1'b0;
            mm_rd_addr   <= '0;
            rp_wr_en     <= 1'b0;
            rp_wr_row    <= '0;
            rp_wr_data   <= '0;
            en_wr_en     <= 1'b0;
            en_wr_addr   <= '0;
            en_wr_data   <= '0;
        end else begin
            // Default pulses low
            mm_rd_en <= 1'b0;
            rp_wr_en <= 1'b0;
            en_wr_en <= 1'b0;

            unique case (state)

                // ── Wait for start ────────────────────────────────────────────
                ST_IDLE: begin
                    if (start) begin
                        clear_idx <= '0;
                        done      <= 1'b0;
                        state     <= ST_CLEAR;
                    end
                end

                // ── Clear degree[]/fill[] for all active rows ─────────────────
                ST_CLEAR: begin
                    degree[clear_idx] <= '0;
                    fill[clear_idx]   <= '0;
                    if (clear_idx == COL_IDX_W_P'(active_n - 1)) begin
                        edge_idx <= '0;
                        state    <= ST_P1_READ;
                    end else begin
                        clear_idx <= clear_idx + 1;
                    end
                end

                // ════════════════════════════════════════════════════════════
                // PASS 1 — count degree[i] and degree[j] for every edge
                // ════════════════════════════════════════════════════════════
                ST_P1_READ: begin
                    if (edge_idx == num_edges) begin
                        // All edges counted -> begin prefix sum
                        clear_idx   <= '0;
                        running_sum <= '0;
                        state       <= ST_FIN_INIT;
                    end else begin
                        mm_rd_en   <= 1'b1;
                        mm_rd_addr <= MAIN_ADDR_WP'(MAIN_EDGE_BASE + edge_idx);
                        state      <= ST_P1_WAIT;
                    end
                end

                ST_P1_WAIT: state <= ST_P1_ACC;

                ST_P1_ACC: begin
                    if (mm_rd_valid) begin
                        // NOTE: assumes no self-loops (r != c), which holds for
                        // all standard Gset instances. If r==c the two
                        // non-blocking writes below would both target
                        // degree[r], and only the LAST one in program order
                        // committing is not guaranteed across simulators --
                        // this is intentionally left unguarded since Gset
                        // graphs never contain self-loops.
                        automatic logic [COL_IDX_W_P-1:0] r, c;
                        r = unpack_row(mm_rd_data);
                        c = unpack_col(mm_rd_data);
                        degree[r] <= degree[r] + 1;
                        degree[c] <= degree[c] + 1;
                    end
                    edge_idx <= edge_idx + 1;
                    state    <= ST_P1_READ;
                end

                // ════════════════════════════════════════════════════════════
                // FINALIZE — prefix sum: rstart[i] = sum(degree[0..i-1])
                // ════════════════════════════════════════════════════════════
                ST_FIN_INIT: state <= ST_FIN_STEP;

                ST_FIN_STEP: begin
                    rstart[clear_idx] <= running_sum;
                    running_sum       <= running_sum + CSR_ADDR_W_P'(degree[clear_idx]);
                    state             <= ST_FIN_WRROWPTR;
                end

                ST_FIN_WRROWPTR: begin
                    // Write row_ptr[clear_idx] = {rstart[clear_idx], degree[clear_idx]}
                    rp_wr_en   <= 1'b1;
                    rp_wr_row  <= clear_idx;
                    rp_wr_data <= {rstart[clear_idx], ROW_COUNT_WP'(degree[clear_idx])};

                    if (clear_idx == COL_IDX_W_P'(active_n - 1)) begin
                        // Prefix sum complete; running_sum now = total nnz
                        nnz_total <= running_sum;
                        edge_idx  <= '0;
                        state     <= ST_P2_READ;
                    end else begin
                        clear_idx <= clear_idx + 1;
                        state     <= ST_FIN_STEP;
                    end
                end

                // ════════════════════════════════════════════════════════════
                // PASS 2 — place each edge's two directed entries
                // ════════════════════════════════════════════════════════════
                ST_P2_READ: begin
                    if (edge_idx == num_edges) begin
                        state <= ST_DONE;
                    end else begin
                        mm_rd_en   <= 1'b1;
                        mm_rd_addr <= MAIN_ADDR_WP'(MAIN_EDGE_BASE + edge_idx);
                        state      <= ST_P2_WAIT;
                    end
                end

                ST_P2_WAIT: state <= ST_P2_PLACE_FWD;

                // Place forward entry: row i -> col j
                ST_P2_PLACE_FWD: begin
                    if (mm_rd_valid) begin
                        automatic logic [COL_IDX_W_P-1:0] row_v, col_v;
                        automatic logic [CSR_ADDR_W_P-1:0] dest_addr;
                        row_v = unpack_row(mm_rd_data);
                        col_v = unpack_col(mm_rd_data);

                        e_row  <= row_v;
                        e_col  <= col_v;
                        e_addr <= MAIN_ADDR_WP'(MAIN_EDGE_BASE + edge_idx);

                        dest_addr  = CSR_ADDR_W_P'(rstart[row_v] + fill[row_v]);
                        en_wr_en   <= 1'b1;
                        en_wr_addr <= dest_addr;
                        en_wr_data <= {col_v, MAIN_ADDR_WP'(MAIN_EDGE_BASE + edge_idx)};
                        fill[row_v] <= fill[row_v] + 1;
                    end
                    state <= ST_P2_PLACE_REV;
                end

                // Place reverse entry: row j -> col i  (same edge, symmetric)
                ST_P2_PLACE_REV: begin
                    en_wr_en   <= 1'b1;
                    en_wr_addr <= CSR_ADDR_W_P'(rstart[e_col] + fill[e_col]);
                    en_wr_data <= {e_row, e_addr};
                    fill[e_col] <= fill[e_col] + 1;

                    edge_idx <= edge_idx + 1;
                    state    <= ST_P2_READ;
                end

                // ── Done ──────────────────────────────────────────────────────
                ST_DONE: begin
                    done <= 1'b1;
                    if (start) begin
                        // Allow re-scan for a new problem instance
                        clear_idx <= '0;
                        done      <= 1'b0;
                        state     <= ST_CLEAR;
                    end
                end

                default: state <= ST_IDLE;

            endcase
        end
    end

endmodule

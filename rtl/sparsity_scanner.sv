// =============================================================================
// sparsity_scanner.sv
//
// One-shot FSM that runs at startup (before the dSB algorithm begins).
// It walks every element of the J matrix in main_memory, identifies non-zero
// entries, and builds the CSR index store (row_ptr + entry arrays).
//
// Operation sequence
// ------------------
//  1. Assert start=1 for one cycle after main_memory has been fully loaded.
//  2. The scanner issues N*N reads to main_memory (one per J[i][j] word).
//  3. For each non-zero J[i][j]:
//       - Writes a CSR entry: {col=j, main_addr = i*N+j} to entry_mem[nnz_ptr]
//       - Increments nnz_ptr
//  4. When column j==N-1 for row i, writes row_ptr[i] = {row_start, row_count}
//  5. After all N rows, asserts done=1 and holds it.
//
// Timing
// ------
//  main_memory has 1-cycle read latency (registered output).
//  State sequence per element:
//    SCAN_READ   → issue rd_en + rd_addr
//    SCAN_WAIT   → absorb the 1-cycle BRAM latency
//    SCAN_EVAL   → rd_valid=1; check for non-zero; optionally write CSR entry
//    SCAN_NEXTCOL→ advance column (or go to SCAN_NEXTROW if j==N-1)
//    SCAN_NEXTROW→ write row_ptr, advance row, back to SCAN_READ
//
// Outputs
// -------
//  To main_memory   : rd_en, rd_addr
//  To csr_index_store:
//    rp_wr_en, rp_wr_row, rp_wr_data   (row pointer writes)
//    en_wr_en, en_wr_addr, en_wr_data  (entry writes)
//  done  : high after the full matrix has been scanned
//  nnz_total : total number of non-zero entries found
// =============================================================================

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
    input  logic start,          // pulse high for one cycle to begin scan

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
    output logic [CSR_ADDR_W_P-1:0]  nnz_total    // total non-zeros found
);

    // ── Internal registers ────────────────────────────────────────────────────
    scan_state_t state;

    logic [COL_IDX_W_P-1:0]  row_i;       // current row   (0..N-1)
    logic [COL_IDX_W_P-1:0]  col_j;       // current col   (0..N-1)
    logic [CSR_ADDR_W_P-1:0] nnz_ptr;     // flat write pointer into entry_mem
    logic [CSR_ADDR_W_P-1:0] row_start;   // entry_mem address where this row began

    // ── Combinational: main memory address of J[row_i][col_j] ─────────────────
    logic [MAIN_ADDR_WP-1:0] j_addr;
    assign j_addr = MAIN_ADDR_WP'(row_i * N_P + col_j);

    // ── Combinational: non-zero test ──────────────────────────────────────────
    // J values are IC_BITS_P-bit signed words stored in MAIN_WORD_WP-wide cells.
    // Non-zero iff the lower IC_BITS_P bits are not all zero.
    logic is_nonzero;
    assign is_nonzero = |mm_rd_data[IC_BITS_P-1:0];

    // ── FSM ───────────────────────────────────────────────────────────────────
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state      <= SCAN_IDLE;
            row_i      <= '0;
            col_j      <= '0;
            nnz_ptr    <= '0;
            row_start  <= '0;
            nnz_total  <= '0;
            done       <= 1'b0;

            mm_rd_en   <= 1'b0;
            mm_rd_addr <= '0;

            rp_wr_en   <= 1'b0;
            rp_wr_row  <= '0;
            rp_wr_data <= '0;

            en_wr_en   <= 1'b0;
            en_wr_addr <= '0;
            en_wr_data <= '0;

        end else begin
            // Default pulse signals low every cycle
            mm_rd_en  <= 1'b0;
            rp_wr_en  <= 1'b0;
            en_wr_en  <= 1'b0;

            unique case (state)

                // ── Wait for start pulse ──────────────────────────────────────
                SCAN_IDLE: begin
                    if (start) begin
                        row_i     <= '0;
                        col_j     <= '0;
                        nnz_ptr   <= '0;
                        row_start <= '0;
                        done      <= 1'b0;
                        state     <= SCAN_READ;
                    end
                end

                // ── Issue read to main_memory ─────────────────────────────────
                SCAN_READ: begin
                    mm_rd_en   <= 1'b1;
                    mm_rd_addr <= j_addr;
                    state      <= SCAN_WAIT;
                end

                // ── Wait one cycle for BRAM output register ───────────────────
                SCAN_WAIT: begin
                    // mm_rd_valid will be high next cycle
                    state <= SCAN_EVAL;
                end

                // ── Evaluate the returned word ────────────────────────────────
                SCAN_EVAL: begin
                    if (mm_rd_valid && is_nonzero) begin
                        // Write CSR entry: {col_j, j_addr}
                        en_wr_en   <= 1'b1;
                        en_wr_addr <= nnz_ptr;
                        en_wr_data <= {col_j, j_addr};
                        nnz_ptr    <= nnz_ptr + 1;
                    end

                    // Decide next: last column of this row?
                    if (col_j == COL_IDX_W_P'(N_P - 1))
                        state <= SCAN_NEXTROW;
                    else
                        state <= SCAN_NEXTCOL;
                end

                // ── Advance column ────────────────────────────────────────────
                SCAN_NEXTCOL: begin
                    col_j <= col_j + 1;
                    state <= SCAN_READ;
                end

                // ── Close row: write row_ptr, check if last row ───────────────
                SCAN_NEXTROW: begin
                    // row_count for this row = nnz_ptr - row_start
                    // (nnz_ptr may have just been incremented in SCAN_EVAL,
                    //  so we use the post-increment value)
                    rp_wr_en   <= 1'b1;
                    rp_wr_row  <= row_i;
                    rp_wr_data <= {row_start,
                                   ROW_COUNT_WP'(nnz_ptr - row_start)};

                    if (row_i == COL_IDX_W_P'(N_P - 1)) begin
                        // Last row done
                        nnz_total <= nnz_ptr;
                        state     <= SCAN_DONE;
                    end else begin
                        // Advance to next row
                        row_start <= nnz_ptr;
                        row_i     <= row_i + 1;
                        col_j     <= '0;
                        state     <= SCAN_READ;
                    end
                end

                // ── Done ──────────────────────────────────────────────────────
                SCAN_DONE: begin
                    done <= 1'b1;
                    // Remain here until reset or a new start pulse
                    if (start) begin
                        // Allow re-scan (e.g. if problem is reloaded)
                        row_i     <= '0;
                        col_j     <= '0;
                        nnz_ptr   <= '0;
                        row_start <= '0;
                        done      <= 1'b0;
                        state     <= SCAN_READ;
                    end
                end

                default: state <= SCAN_IDLE;

            endcase
        end
    end

endmodule

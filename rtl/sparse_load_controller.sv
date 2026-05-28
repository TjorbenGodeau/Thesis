// =============================================================================
// sparse_load_controller.sv
//
// Central orchestration FSM for one oscillator update in the sparse dSB machine.
//
// For each oscillator i (one call per dSB step), this controller:
//   1. Reads row_ptr[i] from csr_index_store → gets {start_addr, count}
//   2. Loops over the row's non-zeros in K_MAX-sized chunks:
//      a. For each entry in the chunk:
//           - Read CSR entry[a] → {col_j, main_addr}
//           - Read J value from main_memory[main_addr]
//           - Write {col_j, J, sign(x_j)} into sparse_compute_memory slot k
//      b. Precharge the 8T tile
//      c. Drive rwl[k] = sign(x_i) for all filled slots
//         (since we're computing J[i][*] · sign(x_*), the RWL is sign(x_i)
//          and the stored value is sign(x_j) — XNOR fires on equality)
//      d. Trigger dotprod_phase2_sparse (accumulate=0 first chunk, 1 after)
//      e. Wait for dotprod done
//   3. After all chunks, asserts done=1.
//      Caller can then read Jx_i from dotprod_phase2_sparse.
//
// Signals requiring attention
// ---------------------------
//  osc_signs [N-1:0]   : live sign vector of ALL oscillators.
//                        Used to pick sign(x_j) when loading a slot.
//                        Must be held stable for the duration of one call.
//
//  sign_xi             : sign of oscillator i (the target).
//                        = osc_signs[osc_idx].  Wired externally or derived.
//
//  rwl [K_MAX-1:0]     : per-slot RWL driven to sparse_compute_memory.
//                        During compute: all slots get sign_xi.
//                        Unused slots get 0 (no discharge).
//
// Chunking example (N=10, K_MAX=4, row i has 7 non-zeros)
// --------------------------------------------------------
//  Chunk 0: slots 0..3 filled, valid_count=4, accumulate=0, last_chunk=0
//  Chunk 1: slots 0..2 filled, valid_count=3, accumulate=1, last_chunk=1
//            (slot 3 left empty / don't-care, rwl[3]=0)
//
// Latency per entry (inside a chunk)
// -----------------------------------
//  SLC_LOAD_ENTRY : issue CSR read
//  SLC_WAIT_ENTRY : absorb 1-cycle CSR BRAM latency
//  SLC_READ_J     : issue main_memory read (using main_addr from CSR entry)
//  SLC_WAIT_J     : absorb 1-cycle main_memory BRAM latency
//  SLC_WRITE_TILE : write col_idx + J into sparse_compute_memory
//  SLC_NEXT_ENTRY : advance slot/entry pointer
//  (6 cycles per entry)
//
// After all K_MAX (or fewer) entries are loaded:
//  SLC_TRIGGER_PH2: assert dotprod start (1 cycle)
//  SLC_WAIT_PH2   : wait for dotprod done (1 cycle in current dotprod impl)
// =============================================================================

module sparse_load_controller
    import dsb_pkg::*;
    import dsb_sparse_pkg::*;
#(
    parameter int unsigned N_P          = N,
    parameter int unsigned K_MAX_P      = K_MAX,
    parameter int unsigned IC_BITS_P    = IC_BITS,
    parameter int unsigned COL_IDX_W_P  = COL_IDX_W,
    parameter int unsigned MAIN_ADDR_WP = MAIN_ADDR_W,
    parameter int unsigned MAIN_WORD_WP = MAIN_WORD_W,
    parameter int unsigned CSR_ADDR_W_P = CSR_ADDR_W,
    parameter int unsigned CSR_ENTRY_WP = CSR_ENTRY_W,
    parameter int unsigned ROW_PTR_W_P  = ROW_PTR_W,
    parameter int unsigned ROW_COUNT_WP = ROW_COUNT_W,
    parameter int unsigned K_CNT_W_P    = $clog2(K_MAX + 1),
    parameter int unsigned XY_W_P       = XY_W
)(
    input  logic clk,
    input  logic rst_n,

    // ── Trigger ───────────────────────────────────────────────────────────────
    input  logic                        start,       // pulse to process osc_idx
    input  logic [COL_IDX_W_P-1:0]      osc_idx,     // which oscillator i to process
    output logic                        done,        // row fully processed

    // ── Live oscillator sign vector ───────────────────────────────────────────
    input  logic [N_P-1:0]              osc_signs,   // sign(x_j) for all j

    // ── csr_index_store read ports ────────────────────────────────────────────
    output logic                        rp_rd_en,
    output logic [COL_IDX_W_P-1:0]      rp_rd_row,
    input  logic [ROW_PTR_W_P-1:0]      rp_rd_data,
    input  logic                        rp_rd_valid,

    output logic                        en_rd_en,
    output logic [CSR_ADDR_W_P-1:0]     en_rd_addr,
    input  logic [CSR_ENTRY_WP-1:0]     en_rd_data,
    input  logic                        en_rd_valid,

    // ── main_memory read port ─────────────────────────────────────────────────
    output logic                        mm_rd_en,
    output logic [MAIN_ADDR_WP-1:0]     mm_rd_addr,
    input  logic [MAIN_WORD_WP-1:0]     mm_rd_data,
    input  logic                        mm_rd_valid,

    // ── sparse_compute_memory write port ─────────────────────────────────────
    output logic                        scm_wr_en,
    output logic [$clog2(K_MAX)-1:0]    scm_wr_slot,
    output logic [COL_IDX_W_P-1:0]      scm_wr_col_idx,
    output logic [IC_BITS_P-1:0]        scm_wr_J,
    output logic                        scm_wr_sign_j,

    // ── sparse_compute_memory compute port ────────────────────────────────────
    output logic                        scm_precharge,
    output logic [K_MAX_P-1:0]          scm_rwl,      // one per slot

    // ── dotprod_phase2_sparse control ─────────────────────────────────────────
    output logic                        ph2_start,
    output logic                        ph2_accumulate,
    output logic                        ph2_last_chunk,
    output logic [K_CNT_W_P-1:0]        ph2_valid_count,
    output logic                        ph2_sign_xi,
    input  logic                        ph2_done
);

    // ── Internal registers ────────────────────────────────────────────────────
    slc_state_t state;

    logic [COL_IDX_W_P-1:0]   cur_osc;          // latched osc_idx
    logic [CSR_ADDR_W_P-1:0]  csr_start;         // start of this row in entry_mem
    logic [ROW_COUNT_WP-1:0]  row_count;          // total non-zeros in this row
    logic [ROW_COUNT_WP-1:0]  entries_remaining;  // entries left to process
    logic [CSR_ADDR_W_P-1:0]  entry_ptr;          // current position in entry_mem
    logic [$clog2(K_MAX)-1:0] slot_ptr;           // current slot in compute memory
    logic [K_CNT_W_P-1:0]     chunk_valid;        // valid entries filled this chunk
    logic                      first_chunk;        // 0 after first chunk

    // Latched values from memory reads
    logic [COL_IDX_W_P-1:0]   latch_col_j;
    logic [MAIN_ADDR_WP-1:0]  latch_main_addr;
    logic [IC_BITS_P-1:0]     latch_J;

    // ── Combinational: sign of current oscillator i ───────────────────────────
    logic sign_xi;
    assign sign_xi = osc_signs[cur_osc];

    // ── Combinational: is this the last chunk? ────────────────────────────────
    // Last chunk when entries_remaining <= K_MAX_P
    logic is_last_chunk;
    assign is_last_chunk = (entries_remaining <= ROW_COUNT_WP'(K_MAX_P));

    // ── FSM ───────────────────────────────────────────────────────────────────
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state             <= SLC_IDLE;
            cur_osc           <= '0;
            csr_start         <= '0;
            row_count         <= '0;
            entries_remaining <= '0;
            entry_ptr         <= '0;
            slot_ptr          <= '0;
            chunk_valid       <= '0;
            first_chunk       <= 1'b1;
            latch_col_j       <= '0;
            latch_main_addr   <= '0;
            latch_J           <= '0;
            done              <= 1'b0;

            rp_rd_en   <= 1'b0;
            rp_rd_row  <= '0;
            en_rd_en   <= 1'b0;
            en_rd_addr <= '0;
            mm_rd_en   <= 1'b0;
            mm_rd_addr <= '0;
            scm_wr_en  <= 1'b0;
            scm_precharge <= 1'b0;
            scm_rwl    <= '0;
            ph2_start  <= 1'b0;

        end else begin
            // Default: pulses low
            rp_rd_en      <= 1'b0;
            en_rd_en      <= 1'b0;
            mm_rd_en      <= 1'b0;
            scm_wr_en     <= 1'b0;
            scm_precharge <= 1'b0;
            ph2_start     <= 1'b0;
            done          <= 1'b0;

            unique case (state)

                // ── Wait for start ────────────────────────────────────────────
                SLC_IDLE: begin
                    if (start) begin
                        cur_osc     <= osc_idx;
                        first_chunk <= 1'b1;
                        done        <= 1'b0;
                        // Read row_ptr[osc_idx]
                        rp_rd_en  <= 1'b1;
                        rp_rd_row <= osc_idx;
                        state     <= SLC_LOAD_PTR;
                    end
                end

                // ── Wait for row_ptr BRAM latency ─────────────────────────────
                SLC_LOAD_PTR: begin
                    state <= SLC_WAIT_PTR;
                end

                SLC_WAIT_PTR: begin
                    if (rp_rd_valid) begin
                        // Unpack {start_addr, count}
                        csr_start  <= rp_rd_data[ROW_PTR_W_P-1 : ROW_COUNT_WP];
                        row_count  <= rp_rd_data[ROW_COUNT_WP-1 : 0];
                        entry_ptr  <= rp_rd_data[ROW_PTR_W_P-1 : ROW_COUNT_WP];
                        entries_remaining <= rp_rd_data[ROW_COUNT_WP-1 : 0];
                        slot_ptr   <= '0;
                        chunk_valid<= '0;
                        state      <= SLC_CHUNK_START;
                    end
                end

                // ── Begin a new chunk ─────────────────────────────────────────
                SLC_CHUNK_START: begin
                    slot_ptr    <= '0;
                    chunk_valid <= '0;
                    if (entries_remaining == '0) begin
                        // Row had zero non-zeros — nothing to compute
                        done  <= 1'b1;
                        state <= SLC_DONE;
                    end else begin
                        // Issue first CSR entry read of this chunk
                        en_rd_en   <= 1'b1;
                        en_rd_addr <= entry_ptr;
                        state      <= SLC_LOAD_ENTRY;
                    end
                end

                // ── Wait for CSR entry BRAM latency ──────────────────────────
                SLC_LOAD_ENTRY: begin
                    state <= SLC_WAIT_ENTRY;
                end

                SLC_WAIT_ENTRY: begin
                    if (en_rd_valid) begin
                        // Unpack {col_j, main_addr}
                        latch_col_j     <= en_rd_data[CSR_ENTRY_WP-1 : MAIN_ADDR_WP];
                        latch_main_addr <= en_rd_data[MAIN_ADDR_WP-1 : 0];
                        // Issue main_memory read for J value
                        mm_rd_en   <= 1'b1;
                        mm_rd_addr <= en_rd_data[MAIN_ADDR_WP-1 : 0];
                        state      <= SLC_READ_J;
                    end
                end

                // ── Wait for main_memory BRAM latency ────────────────────────
                SLC_READ_J: begin
                    state <= SLC_WAIT_J;
                end

                SLC_WAIT_J: begin
                    if (mm_rd_valid) begin
                        latch_J <= mm_rd_data[IC_BITS_P-1:0];
                        state   <= SLC_WRITE_TILE;
                    end
                end

                // ── Write entry into sparse_compute_memory ────────────────────
                SLC_WRITE_TILE: begin
                    scm_wr_en     <= 1'b1;
                    scm_wr_slot   <= slot_ptr;
                    scm_wr_col_idx<= latch_col_j;
                    scm_wr_J      <= latch_J;
                    scm_wr_sign_j <= osc_signs[latch_col_j]; // live sign(x_j)

                    entries_remaining <= entries_remaining - 1;
                    entry_ptr         <= entry_ptr + 1;
                    chunk_valid       <= chunk_valid + 1;

                    state <= SLC_NEXT_ENTRY;
                end

                // ── Decide: more slots in this chunk, or trigger dotprod? ─────
                SLC_NEXT_ENTRY: begin
                    // Chunk full (slot_ptr == K_MAX_P-1), OR last entry of row
                    if (slot_ptr == $clog2(K_MAX)'(K_MAX_P - 1) ||
                        entries_remaining == '0) begin
                        // Trigger dotprod for this chunk — precharge first
                        scm_precharge <= 1'b1;
                        for (int k = 0; k < K_MAX_P; k++) begin
                            scm_rwl[k] <= (k <= int'(slot_ptr)) ? sign_xi : 1'b0;
                        end
                        state <= SLC_TRIGGER_PH2;
                    end else begin
                        // More entries to load: advance slot and fetch next entry
                        slot_ptr   <= slot_ptr + 1;
                        en_rd_en   <= 1'b1;
                        en_rd_addr <= entry_ptr;
                        state      <= SLC_LOAD_ENTRY;
                    end
                end

                // ── Trigger dotprod_phase2_sparse ─────────────────────────────
                SLC_TRIGGER_PH2: begin
                    ph2_start       <= 1'b1;
                    ph2_accumulate  <= ~first_chunk;
                    ph2_last_chunk  <= is_last_chunk;
                    // On the last chunk supply total row_count so that the
                    // sign-correction in dotprod_phase2_sparse accounts for
                    // ALL non-zero terms (not just those in the final chunk).
                    ph2_valid_count <= is_last_chunk
                                       ? K_CNT_W_P'(row_count)
                                       : K_CNT_W_P'(chunk_valid);
                    ph2_sign_xi     <= sign_xi;
                    first_chunk     <= 1'b0;
                    state           <= SLC_WAIT_PH2;
                end

                // ── Wait for dotprod to complete ──────────────────────────────
                SLC_WAIT_PH2: begin
                    if (ph2_done) begin
                        state <= SLC_NEXT_CHUNK;
                    end
                end

                // ── Move to next chunk or finish ──────────────────────────────
                SLC_NEXT_CHUNK: begin
                    if (entries_remaining == '0) begin
                        done  <= 1'b1;
                        state <= SLC_DONE;
                    end else begin
                        slot_ptr    <= '0;
                        chunk_valid <= '0;
                        // Issue first CSR entry read of next chunk
                        en_rd_en   <= 1'b1;
                        en_rd_addr <= entry_ptr;
                        state      <= SLC_LOAD_ENTRY;
                    end
                end

                // ── Done ──────────────────────────────────────────────────────
                SLC_DONE: begin
                    // Hold done high for one cycle then return to idle
                    state <= SLC_IDLE;
                end

                default: state <= SLC_IDLE;

            endcase
        end
    end

endmodule

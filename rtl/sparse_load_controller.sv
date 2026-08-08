// =============================================================================
// sparse_load_controller.sv (fixed valid_count)
// =============================================================================
`timescale 1ns/1ps
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
    parameter int unsigned K_CNT_W_P    = $clog2(K_MAX_P + 1),
    parameter int unsigned XY_W_P       = XY_W
)(
    input  logic clk,
    input  logic rst_n,
    input  logic                        start,
    input  logic [COL_IDX_W_P-1:0]      osc_idx,
    output logic                        done,
    input  logic [N_P-1:0]              osc_signs,
    output logic                        rp_rd_en,
    output logic [COL_IDX_W_P-1:0]      rp_rd_row,
    input  logic [ROW_PTR_W_P-1:0]      rp_rd_data,
    input  logic                        rp_rd_valid,
    output logic                        en_rd_en,
    output logic [CSR_ADDR_W_P-1:0]     en_rd_addr,
    input  logic [CSR_ENTRY_WP-1:0]     en_rd_data,
    input  logic                        en_rd_valid,
    output logic                        mm_rd_en,
    output logic [MAIN_ADDR_WP-1:0]     mm_rd_addr,
    input  logic [MAIN_WORD_WP-1:0]     mm_rd_data,
    input  logic                        mm_rd_valid,
    output logic                        scm_wr_en,
    output logic [$clog2(K_MAX_P)-1:0]  scm_wr_slot,
    output logic [COL_IDX_W_P-1:0]      scm_wr_col_idx,
    output logic [IC_BITS_P-1:0]        scm_wr_J,
    output logic                        scm_wr_sign_j,
    output logic                        scm_precharge,
    output logic [K_MAX_P-1:0]          scm_rwl,
    output logic                        ph2_start,
    output logic                        ph2_accumulate,
    output logic                        ph2_last_chunk,
    output logic [K_CNT_W_P-1:0]        ph2_valid_count,
    output logic                        ph2_sign_xi,
    input  logic                        ph2_done
);

    typedef enum logic [4:0] {
        SLC_IDLE,
        SLC_LOAD_PTR,
        SLC_WAIT_PTR,
        SLC_CHUNK_START,
        SLC_LOAD_ENTRY,
        SLC_WAIT_ENTRY,
        SLC_READ_J,
        SLC_WAIT_J,
        SLC_WRITE_TILE,
        SLC_NEXT_ENTRY,
        SLC_PRECLEAR,
        SLC_TRIGGER_PH2,
        SLC_WAIT_PH2,
        SLC_NEXT_CHUNK,
        SLC_DONE
    } slc_state_t;

    slc_state_t state;
    logic [COL_IDX_W_P-1:0]   cur_osc;
    logic [CSR_ADDR_W_P-1:0]  csr_start;
    logic [ROW_COUNT_WP-1:0]  row_count;
    logic [ROW_COUNT_WP-1:0]  entries_remaining;
    logic [CSR_ADDR_W_P-1:0]  entry_ptr;
    logic [$clog2(K_MAX_P)-1:0] slot_ptr;
    logic [K_CNT_W_P-1:0]     chunk_valid;
    logic                      first_chunk;
    logic [COL_IDX_W_P-1:0]   latch_col_j;
    logic [MAIN_ADDR_WP-1:0]  latch_main_addr;
    logic [IC_BITS_P-1:0]     latch_J;

    logic sign_xi;
    assign sign_xi = osc_signs[cur_osc];

    logic is_last_chunk;
    assign is_last_chunk = (entries_remaining <= ROW_COUNT_WP'(K_MAX_P));

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
            rp_rd_en      <= 1'b0;
            en_rd_en      <= 1'b0;
            mm_rd_en      <= 1'b0;
            scm_wr_en     <= 1'b0;
            scm_precharge <= 1'b0;
            ph2_start     <= 1'b0;
            done          <= 1'b0;

            unique case (state)
                SLC_IDLE: begin
                    if (start) begin
                        cur_osc     <= osc_idx;
                        first_chunk <= 1'b1;
                        done        <= 1'b0;
                        rp_rd_en  <= 1'b1;
                        rp_rd_row <= osc_idx;
                        state     <= SLC_LOAD_PTR;
                    end
                end
                SLC_LOAD_PTR: begin
                    state <= SLC_WAIT_PTR;
                end
                SLC_WAIT_PTR: begin
                    if (rp_rd_valid) begin
                        csr_start  <= rp_rd_data[ROW_PTR_W_P-1 : ROW_COUNT_WP];
                        row_count  <= rp_rd_data[ROW_COUNT_WP-1 : 0];
                        entry_ptr  <= rp_rd_data[ROW_PTR_W_P-1 : ROW_COUNT_WP];
                        entries_remaining <= rp_rd_data[ROW_COUNT_WP-1 : 0];
                        slot_ptr   <= '0;
                        chunk_valid<= '0;
                        state      <= SLC_CHUNK_START;
                    end
                end
                SLC_CHUNK_START: begin
                    slot_ptr    <= '0;
                    chunk_valid <= '0;
                    if (entries_remaining == '0) begin
                        done  <= 1'b1;
                        state <= SLC_DONE;
                    end else begin
                        en_rd_en   <= 1'b1;
                        en_rd_addr <= entry_ptr;
                        state      <= SLC_LOAD_ENTRY;
                    end
                end
                SLC_LOAD_ENTRY: begin
                    state <= SLC_WAIT_ENTRY;
                end
                SLC_WAIT_ENTRY: begin
                    if (en_rd_valid) begin
                        latch_col_j     <= en_rd_data[CSR_ENTRY_WP-1 : MAIN_ADDR_WP];
                        latch_main_addr <= en_rd_data[MAIN_ADDR_WP-1 : 0];
                        mm_rd_en   <= 1'b1;
                        mm_rd_addr <= en_rd_data[MAIN_ADDR_WP-1 : 0];
                        state      <= SLC_READ_J;
                    end
                end
                SLC_READ_J: begin
                    state <= SLC_WAIT_J;
                end
                SLC_WAIT_J: begin
                    if (mm_rd_valid) begin
                        latch_J <= mm_rd_data[IC_BITS_P-1:0];
                        state   <= SLC_WRITE_TILE;
                    end
                end
                SLC_WRITE_TILE: begin
                    scm_wr_en     <= 1'b1;
                    scm_wr_slot   <= slot_ptr;
                    scm_wr_col_idx<= latch_col_j;
                    scm_wr_J      <= latch_J;
                    scm_wr_sign_j <= osc_signs[latch_col_j];
                    entries_remaining <= entries_remaining - 1;
                    entry_ptr         <= entry_ptr + 1;
                    chunk_valid       <= chunk_valid + 1;
                    state <= SLC_NEXT_ENTRY;
                end
                SLC_NEXT_ENTRY: begin
                    if (slot_ptr == $clog2(K_MAX_P)'(K_MAX_P - 1) ||
                        entries_remaining == '0) begin
                        scm_precharge <= 1'b1;
                        for (int k = 0; k < K_MAX_P; k++) begin
                            scm_rwl[k] <= (k <= int'(slot_ptr)) ? sign_xi : 1'b0;
                        end
                        state <= SLC_PRECLEAR;
                    end else begin
                        slot_ptr   <= slot_ptr + 1;
                        en_rd_en   <= 1'b1;
                        en_rd_addr <= entry_ptr;
                        state      <= SLC_LOAD_ENTRY;
                    end
                end
                SLC_PRECLEAR: begin
                    state <= SLC_TRIGGER_PH2;
                end
                SLC_TRIGGER_PH2: begin
                    ph2_start       <= 1'b1;
                    ph2_accumulate  <= ~first_chunk;
                    ph2_last_chunk  <= is_last_chunk;
                    ph2_valid_count <= K_CNT_W_P'(chunk_valid);   // FIX: always use chunk size
                    ph2_sign_xi     <= sign_xi;
                    first_chunk     <= 1'b0;
                    state           <= SLC_WAIT_PH2;
                end
                SLC_WAIT_PH2: begin
                    if (ph2_done) begin
                        state <= SLC_NEXT_CHUNK;
                    end
                end
                SLC_NEXT_CHUNK: begin
                    if (entries_remaining == '0) begin
                        done  <= 1'b1;
                        state <= SLC_DONE;
                    end else begin
                        slot_ptr    <= '0;
                        chunk_valid <= '0;
                        en_rd_en   <= 1'b1;
                        en_rd_addr <= entry_ptr;
                        state      <= SLC_LOAD_ENTRY;
                    end
                end
                SLC_DONE: begin
                    state <= SLC_IDLE;
                end
                default: state <= SLC_IDLE;
            endcase
        end
    end
endmodule
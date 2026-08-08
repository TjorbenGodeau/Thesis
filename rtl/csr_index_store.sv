// =============================================================================
// csr_index_store.sv
//
// Compressed Sparse Row (CSR) index store for the sparse dSB Ising machine.
//
// Two independent sub-memories, both with 1-cycle read latency:
//
//  1. row_ptr_mem  [N entries]
//     Address  : row index i  (COL_IDX_W bits, reused as row addr)
//     Data     : {start_addr[CSR_ADDR_W-1:0], count[ROW_COUNT_W-1:0]}
//                = ROW_PTR_W bits total
//     Meaning  : row i's non-zero entries start at CSR entry address
//                'start_addr' and there are 'count' of them.
//
//  2. entry_mem    [MAX_NNZ entries]
//     Address  : flat CSR index  a  (CSR_ADDR_W bits)
//     Data     : {col_idx[COL_IDX_W-1:0], main_addr[MAIN_ADDR_W-1:0]}
//                = CSR_ENTRY_W bits total
//     Meaning  : the non-zero at flat index a is J[row][col_idx], and
//                its value lives at main_addr in main_memory.
//
// Write ports
// -----------
//  Row pointer write : rp_wr_en, rp_wr_row, rp_wr_data
//  Entry write       : en_wr_en, en_wr_addr, en_wr_data
//  (Both are written once by sparsity_scanner at startup.)
//
// Read ports
// ----------
//  Row pointer read  : rp_rd_en, rp_rd_row → rp_rd_data / rp_rd_valid
//  Entry read        : en_rd_en, en_rd_addr → en_rd_data / en_rd_valid
//  (Read repeatedly by sparse_load_controller during the algorithm run.)
//
// Bit-field helpers (pack / unpack in the caller)
// ------------------------------------------------
//  row_ptr  packing : rp_wr_data = {start_addr, count}
//                     start = rp_rd_data[ROW_PTR_W-1 : ROW_COUNT_W]
//                     count = rp_rd_data[ROW_COUNT_W-1 : 0]
//
//  entry    packing : en_wr_data = {col_idx, main_addr}
//                     col   = en_rd_data[CSR_ENTRY_W-1 : MAIN_ADDR_W]
//                     maddr = en_rd_data[MAIN_ADDR_W-1 : 0]
// =============================================================================

`timescale 1ns/1ps
module csr_index_store
    import dsb_pkg::*;
    import dsb_sparse_pkg::*;
#(
    parameter int unsigned N_P          = N,
    parameter int unsigned MAX_NNZ_P    = MAX_NNZ,
    parameter int unsigned COL_IDX_W_P  = COL_IDX_W,
    parameter int unsigned MAIN_ADDR_WP = MAIN_ADDR_W,
    parameter int unsigned CSR_ADDR_W_P = CSR_ADDR_W,
    parameter int unsigned CSR_ENTRY_WP = CSR_ENTRY_W,
    parameter int unsigned ROW_PTR_W_P  = ROW_PTR_W,
    parameter int unsigned ROW_COUNT_WP = ROW_COUNT_W
)(
    input  logic clk,
    input  logic rst_n,

    // ── Row-pointer write port ────────────────────────────────────────────────
    input  logic                    rp_wr_en,
    input  logic [COL_IDX_W_P-1:0]  rp_wr_row,    // row index i
    input  logic [ROW_PTR_W_P-1:0]  rp_wr_data,   // {start_addr, count}

    // ── Row-pointer read port ─────────────────────────────────────────────────
    input  logic                    rp_rd_en,
    input  logic [COL_IDX_W_P-1:0]  rp_rd_row,
    output logic [ROW_PTR_W_P-1:0]  rp_rd_data,   // valid one cycle later
    output logic                    rp_rd_valid,

    // ── Entry write port ──────────────────────────────────────────────────────
    input  logic                    en_wr_en,
    input  logic [CSR_ADDR_W_P-1:0] en_wr_addr,
    input  logic [CSR_ENTRY_WP-1:0] en_wr_data,   // {col_idx, main_addr}

    // ── Entry read port ───────────────────────────────────────────────────────
    input  logic                    en_rd_en,
    input  logic [CSR_ADDR_W_P-1:0] en_rd_addr,
    output logic [CSR_ENTRY_WP-1:0] en_rd_data,   // valid one cycle later
    output logic                    en_rd_valid
);

    // ── Row-pointer memory ────────────────────────────────────────────────────
    logic [ROW_PTR_W_P-1:0] row_ptr_mem [0:N_P-1];

    always_ff @(posedge clk) begin
        if (rp_wr_en)
            row_ptr_mem[rp_wr_row] <= rp_wr_data;
    end

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            rp_rd_data  <= '0;
            rp_rd_valid <= 1'b0;
        end else begin
            rp_rd_valid <= 1'b0;
            if (rp_rd_en) begin
                rp_rd_data  <= row_ptr_mem[rp_rd_row];
                rp_rd_valid <= 1'b1;
            end
        end
    end

    // ── Entry memory ──────────────────────────────────────────────────────────
    logic [CSR_ENTRY_WP-1:0] entry_mem [0:MAX_NNZ_P-1];

    always_ff @(posedge clk) begin
        if (en_wr_en)
            entry_mem[en_wr_addr] <= en_wr_data;
    end

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            en_rd_data  <= '0;
            en_rd_valid <= 1'b0;
        end else begin
            en_rd_valid <= 1'b0;
            if (en_rd_en) begin
                en_rd_data  <= entry_mem[en_rd_addr];
                en_rd_valid <= 1'b1;
            end
        end
    end

endmodule

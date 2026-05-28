// =============================================================================
// main_memory.sv
//
// Word-addressable main memory for the sparse dSB Ising machine.
//
// Memory map
// ----------
//  [0 .. N*N-1]            J[i][j]  at address i*N + j   (IC_BITS signed)
//  [N*N .. N*N+N-1]        x_init[i]                     (XY_W signed)
//  [N*N+N .. N*N+2N-1]     y_init[i]                     (XY_W signed)
//
// All words are stored at MAIN_WORD_W (= XY_W) bits wide.
// J values are sign-extended on write, sign-truncated on read.
//
// Port summary
// ------------
//  Write port  (WR) : single-cycle, synchronous write.
//  Read port   (RD) : registered output (1-cycle latency) — standard BRAM
//                     inference style.  rd_valid pulses the cycle after rd_en.
//
// Testbench hint
// --------------
//  To preload J:   wr_en=1, wr_addr = i*N+j, wr_data = J[i][j] sign-extended
//  To preload x/y: wr_en=1, wr_addr = MAIN_X_BASE+i / MAIN_Y_BASE+i
//  To read J[i][j]: rd_en=1, rd_addr = i*N+j  → data appears on rd_data one
//                   cycle later with rd_valid=1
// =============================================================================

module main_memory
    import dsb_pkg::*;
    import dsb_sparse_pkg::*;
#(
    parameter int unsigned WORDS_P    = MAIN_MEM_WORDS,
    parameter int unsigned WORD_W_P   = MAIN_WORD_W,
    parameter int unsigned ADDR_W_P   = MAIN_ADDR_W
)(
    input  logic                  clk,
    input  logic                  rst_n,

    // ── Write port ────────────────────────────────────────────────────────────
    input  logic                  wr_en,
    input  logic [ADDR_W_P-1:0]   wr_addr,
    input  logic [WORD_W_P-1:0]   wr_data,    // sign-extend J before driving

    // ── Read port ─────────────────────────────────────────────────────────────
    input  logic                  rd_en,
    input  logic [ADDR_W_P-1:0]   rd_addr,
    output logic [WORD_W_P-1:0]   rd_data,    // valid one cycle after rd_en
    output logic                  rd_valid     // single-cycle pulse
);

    // ── Storage array ─────────────────────────────────────────────────────────
    logic [WORD_W_P-1:0] mem [0:WORDS_P-1];

    // ── Write ─────────────────────────────────────────────────────────────────
    always_ff @(posedge clk) begin
        if (wr_en)
            mem[wr_addr] <= wr_data;
    end

    // ── Read (registered — matches BRAM output-register mode) ─────────────────
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            rd_data  <= '0;
            rd_valid <= 1'b0;
        end else begin
            rd_valid <= 1'b0;   // default: no valid data
            if (rd_en) begin
                rd_data  <= mem[rd_addr];
                rd_valid <= 1'b1;
            end
        end
    end

    // ── Address helper functions (for use in testbenches / other modules) ──────
    // J[i][j]   → address  i*N + j
    // x_init[i] → address  MAIN_X_BASE + i
    // y_init[i] → address  MAIN_Y_BASE + i
    //
    // These are just documentation — callers compute the address themselves.

endmodule

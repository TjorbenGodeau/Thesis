// =============================================================================
// sparse_compute_memory.sv
//
// The compute-side storage for one chunk of K_MAX sparse entries.
// Two physically separate arrays, written together and read/used separately:
//
//  1. col_idx_sram  [K_MAX entries × COL_IDX_W bits]
//     Plain synchronous SRAM.
//     Holds the column index (oscillator index j) for each slot k.
//     Read by sparse_load_controller to determine which sign_j to feed as RWL.
//
//  2. j_bitcell_row [K_MAX entries × IC_BITS rows of dsb_8t_bitcell]
//     8T SRAM compute fabric.
//     Holds the J magnitude bits for each slot k.
//     The RWL for slot k is sign(x_{col_idx[k]}) — driven externally.
//     The XNOR result goes to xnor_phase1 / dotprod_phase2_sparse.
//
// Write interface (from sparse_load_controller)
// ----------------------------------------------
//  wr_en        : write enable (one slot per cycle)
//  wr_slot      : slot index k  (0..K_MAX-1)
//  wr_col_idx   : column index j to store
//  wr_J         : IC_BITS-wide J magnitude bits
//
// Compute interface (driven by core FSM, same as original compute_tile)
// ---------------------------------------------------------------------
//  precharge    : precharge all RBLs
//  rwl[k]       : sign(x_{col_idx[k]}) — one RWL per slot
//                 Caller reads col_idx_sram first, looks up sign of that osc.
//  rbl_J[k*IC_BITS +: IC_BITS] : XNOR output per slot
//  rbl_signs[k] : sign-equality output per slot (sign bitcell)
//
// Read interface for col_idx (combinational, registered output)
// -------------------------------------------------------------
//  rd_slot      : slot to read
//  rd_col_idx   : registered output (1-cycle latency)
//  rd_col_valid : valid strobe
//
// Note: The sign bitcell (one per slot) stores sign(x_j) at write time.
//       rwl[k] = sign(x_i) is the TARGET oscillator's sign (not the stored one).
//       The XNOR fires when stored sign == rwl, i.e. sign(x_j) == sign(x_i).
// =============================================================================

module sparse_compute_memory
    import dsb_pkg::*;
    import dsb_sparse_pkg::*;
#(
    parameter int unsigned K_MAX_P   = K_MAX,
    parameter int unsigned IC_BITS_P = IC_BITS,
    parameter int unsigned COL_W_P   = COL_IDX_W
)(
    input  logic clk,

    // ── Write port (load new chunk) ───────────────────────────────────────────
    input  logic                      wr_en,
    input  logic [$clog2(K_MAX_P)-1:0]  wr_slot,
    input  logic [COL_W_P-1:0]        wr_col_idx,
    input  logic [IC_BITS_P-1:0]      wr_J,
    input  logic                      wr_sign_j,   // sign(x_j) at load time

    // ── Col-index read port ───────────────────────────────────────────────────
    input  logic                      rd_col_en,
    input  logic [$clog2(K_MAX_P)-1:0]  rd_slot,
    output logic [COL_W_P-1:0]        rd_col_idx,
    output logic                      rd_col_valid,

    // ── 8T compute port ───────────────────────────────────────────────────────
    input  logic                      precharge,
    input  logic [K_MAX_P-1:0]        rwl,          // one per slot: sign(x_i)
    output logic [K_MAX_P*IC_BITS_P-1:0] rbl_J,     // XNOR output
    output logic [K_MAX_P-1:0]        rbl_signs      // sign-eq output
);

    // ── Column index SRAM ─────────────────────────────────────────────────────
    logic [COL_W_P-1:0] col_idx_mem [0:K_MAX_P-1];

    always_ff @(posedge clk) begin
        if (wr_en)
            col_idx_mem[wr_slot] <= wr_col_idx;
    end

    always_ff @(posedge clk) begin
        rd_col_valid <= 1'b0;
        if (rd_col_en) begin
            rd_col_idx   <= col_idx_mem[rd_slot];
            rd_col_valid <= 1'b1;
        end
    end

    // ── 8T J bitcell array ────────────────────────────────────────────────────
    // K_MAX_P slots × IC_BITS_P bits each
    // wwl for slot k fires when wr_en && wr_slot==k
    for (genvar k = 0; k < K_MAX_P; k++) begin : gen_slot
        logic wwl_k;
        assign wwl_k = wr_en && (wr_slot == $clog2(K_MAX_P)'(k));

        for (genvar b = 0; b < IC_BITS_P; b++) begin : gen_bit
            dsb_8t_bitcell u_J (
                .clk       (clk),
                .wwl       (wwl_k),
                .wbl       (wr_J[b]),
                .rwl       (rwl[k]),
                .precharge (precharge),
                .rbl_out   (rbl_J[k*IC_BITS_P + b])
            );
        end

        // Sign bitcell: stores sign(x_j), compute fires on sign(x_i)==sign(x_j)
        dsb_8t_bitcell u_sign (
            .clk       (clk),
            .wwl       (wwl_k),
            .wbl       (wr_sign_j),
            .rwl       (rwl[k]),
            .precharge (precharge),
            .rbl_out   (rbl_signs[k])
        );
    end

endmodule

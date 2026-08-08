// =============================================================================
// sparse_core.sv
//
// Drop-in replacement for the original "core" module.
// External interface is IDENTICAL to core — same ports, same done timing —
// so sparse_dsb_array can substitute it without changing its own FSM at all.
//
// What changes inside vs the original core
// -----------------------------------------
//  REMOVED : compute_tile   (dense N-wide 8T tile, loaded by wwl/wbl_J)
//  REMOVED : xnor_phase1    (captures dense RBL row)
//  REMOVED : dotprod_phase2 (accumulates dense N-entry dot product)
//
//  ADDED   : sparse_compute_memory  (K_MAX-slot 8T tile)
//  ADDED   : dotprod_phase2_sparse  (chunked accumulator)
//  ADDED   : sparse_load_controller (reads CSR store + main_memory, fills tile)
//
//  KEPT    : update_unit  (phase 3 — unchanged, identical wiring)
//  KEPT    : schedule     (called from sparse_dsb_array above, a_m passed in)
//
// Port differences from original core
// ------------------------------------
//  REMOVED : wwl, wbl_J, wbl_signs
//            (J data now comes from main_memory via the SLC, not a push bus)
//  ADDED   : all sparse subsystem ports (CSR read, main_memory read, osc_signs)
//
// Timing
// ------
//  Original core: 4 cycles (PRECHARGE → PHASE1 → PHASE2 → PHASE3)
//  Sparse core  : variable — SLC takes 6*nnz_i cycles to load + dotprod (1 cycle
//                 per chunk) + update_unit (1 cycle).  done fires when
//                 update_unit fires done, exactly as before.
//
// FSM states
// ----------
//  SC_IDLE      : waiting for start
//  SC_SLC_RUN   : sparse_load_controller running (loads tile + computes Jx_i)
//  SC_UPDATE    : fire update_unit, wait one cycle for done
// =============================================================================

`timescale 1ns/1ps
module sparse_core
    import dsb_pkg::*;
    import dsb_sparse_pkg::*;
#(
    parameter int unsigned N_P        = N,
    parameter int unsigned IC_BITS_P  = IC_BITS,
    parameter int unsigned XY_W_P     = XY_W,
    parameter int unsigned XY_FRAC_P  = XY_FRAC,
    parameter int unsigned A_BITS_P   = A_BITS,
    parameter int unsigned ACCUM_W_P  = ACCUM_W,
    parameter int unsigned PROD_W_P   = PROD_W,
    parameter int unsigned K_MAX_P    = K_MAX,
    parameter int unsigned COL_IDX_WP = COL_IDX_W,
    parameter int unsigned MAIN_AW    = MAIN_ADDR_W,
    parameter int unsigned MAIN_WW    = MAIN_WORD_W,
    parameter int unsigned CSR_AW     = CSR_ADDR_W,
    parameter int unsigned CSR_EW     = CSR_ENTRY_W,
    parameter int unsigned ROW_PTR_WP = ROW_PTR_W,
    parameter int unsigned ROW_CNT_WP = ROW_COUNT_W
)(
    input  logic                         clk,
    input  logic                         rst_n,
    input  logic                         start,

    // ── Current oscillator state (same as original core) ─────────────────────
    input  logic [COL_IDX_WP-1:0]        osc_idx,      // which oscillator i
    input  logic signed [XY_W_P-1:0]     x_i,
    input  logic signed [XY_W_P-1:0]     y_i,

    // ── Schedule (same as original core) ─────────────────────────────────────
    input  logic [A_BITS_P-1:0]          a_m,
    input  logic [XY_FRAC_P-1:0]         dt_fp,
    input  logic [XY_FRAC_P-1:0]         c0_fp,

    // ── Live sign vector of ALL oscillators ──────────────────────────────────
    input  logic [N_P-1:0]               osc_signs,    // signs_out from array

    // ── CSR index store read ports ────────────────────────────────────────────
    output logic                         rp_rd_en,
    output logic [COL_IDX_WP-1:0]        rp_rd_row,
    input  logic [ROW_PTR_WP-1:0]        rp_rd_data,
    input  logic                         rp_rd_valid,

    output logic                         en_rd_en,
    output logic [CSR_AW-1:0]            en_rd_addr,
    input  logic [CSR_EW-1:0]            en_rd_data,
    input  logic                         en_rd_valid,

    // ── Main memory read port ─────────────────────────────────────────────────
    output logic                         mm_rd_en,
    output logic [MAIN_AW-1:0]           mm_rd_addr,
    input  logic [MAIN_WW-1:0]           mm_rd_data,
    input  logic                         mm_rd_valid,

    // ── Outputs (same as original core) ──────────────────────────────────────
    output logic signed [XY_W_P-1:0]    x_i_new,
    output logic signed [XY_W_P-1:0]    y_i_new,
    output logic                         sign_xi_new,
    output logic                         done
);

    // ── Sign conventions (match original core exactly) ────────────────────────
    assign sign_xi_new = ~x_i_new[XY_W_P-1];

    // =========================================================================
    // sparse_load_controller
    // =========================================================================
    logic                        slc_start;
    logic                        slc_done;

    // SLC → sparse_compute_memory
    logic                        scm_wr_en;
    logic [$clog2(K_MAX)-1:0]    scm_wr_slot;
    logic [COL_IDX_WP-1:0]       scm_wr_col_idx;
    logic [IC_BITS_P-1:0]        scm_wr_J;
    logic                        scm_wr_sign_j;
    logic                        scm_precharge;
    logic [K_MAX_P-1:0]          scm_rwl;

    // SLC → dotprod_phase2_sparse
    logic                        ph2_start;
    logic                        ph2_accumulate;
    logic                        ph2_last_chunk;
    logic [$clog2(K_MAX+1)-1:0]  ph2_valid_count;
    logic                        ph2_sign_xi;
    logic                        ph2_done;

    sparse_load_controller #(
        .N_P         (N_P),
        .K_MAX_P     (K_MAX_P),
        .IC_BITS_P   (IC_BITS_P),
        .COL_IDX_W_P (COL_IDX_WP),
        .MAIN_ADDR_WP(MAIN_AW),
        .MAIN_WORD_WP(MAIN_WW),
        .CSR_ADDR_W_P(CSR_AW),
        .CSR_ENTRY_WP(CSR_EW),
        .ROW_PTR_W_P (ROW_PTR_WP),
        .ROW_COUNT_WP(ROW_CNT_WP)
    ) u_slc (
        .clk            (clk),
        .rst_n          (rst_n),
        .start          (slc_start),
        .osc_idx        (osc_idx),
        .done           (slc_done),
        .osc_signs      (osc_signs),
        .rp_rd_en       (rp_rd_en),
        .rp_rd_row      (rp_rd_row),
        .rp_rd_data     (rp_rd_data),
        .rp_rd_valid    (rp_rd_valid),
        .en_rd_en       (en_rd_en),
        .en_rd_addr     (en_rd_addr),
        .en_rd_data     (en_rd_data),
        .en_rd_valid    (en_rd_valid),
        .mm_rd_en       (mm_rd_en),
        .mm_rd_addr     (mm_rd_addr),
        .mm_rd_data     (mm_rd_data),
        .mm_rd_valid    (mm_rd_valid),
        .scm_wr_en      (scm_wr_en),
        .scm_wr_slot    (scm_wr_slot),
        .scm_wr_col_idx (scm_wr_col_idx),
        .scm_wr_J       (scm_wr_J),
        .scm_wr_sign_j  (scm_wr_sign_j),
        .scm_precharge  (scm_precharge),
        .scm_rwl        (scm_rwl),
        .ph2_start      (ph2_start),
        .ph2_accumulate (ph2_accumulate),
        .ph2_last_chunk (ph2_last_chunk),
        .ph2_valid_count(ph2_valid_count),
        .ph2_sign_xi    (ph2_sign_xi),
        .ph2_done       (ph2_done)
    );

    // =========================================================================
    // sparse_compute_memory
    // =========================================================================
    logic [K_MAX_P*IC_BITS_P-1:0] scm_rbl_J;
    logic [K_MAX_P-1:0]           scm_rbl_signs;

    sparse_compute_memory #(
        .K_MAX_P  (K_MAX_P),
        .IC_BITS_P(IC_BITS_P),
        .COL_W_P  (COL_IDX_WP)
    ) u_scm (
        .clk         (clk),
        .wr_en       (scm_wr_en),
        .wr_slot     (scm_wr_slot),
        .wr_col_idx  (scm_wr_col_idx),
        .wr_J        (scm_wr_J),
        .wr_sign_j   (scm_wr_sign_j),
        .rd_col_en   (1'b0),
        .rd_slot     ('0),
        .rd_col_idx  (),
        .rd_col_valid(),
        .precharge   (scm_precharge),
        .rwl         (scm_rwl),
        .rbl_J       (scm_rbl_J),
        .rbl_signs   (scm_rbl_signs)
    );

    // =========================================================================
    // dotprod_phase2_sparse
    // =========================================================================
    logic signed [ACCUM_W_P-1:0] Jx_i_w;

    dotprod_phase2_sparse #(
        .K_MAX_P   (K_MAX_P),
        .IC_BITS_P (IC_BITS_P),
        .ACCUM_W_P (ACCUM_W_P)
    ) u_ph2 (
        .clk         (clk),
        .start       (ph2_start),
        .done        (ph2_done),
        .accumulate  (ph2_accumulate),
        .last_chunk  (ph2_last_chunk),
        .valid_count (ph2_valid_count),
        .sign_xi     (ph2_sign_xi),
        .xnor_J      (scm_rbl_J),
        .sign_eq     (scm_rbl_signs),
        .Jx_i        (Jx_i_w)
    );

    // =========================================================================
    // update_unit (phase 3 — unchanged from original)
    // =========================================================================
    logic                      ph3_start;
    logic signed [XY_W_P-1:0] x_new_w, y_new_w;
    logic                      ph3_done;
    // Latch x_i and y_i at start (same as original core does at S_IDLE→S_PRECHARGE)
    logic signed [XY_W_P-1:0] x_i_r, y_i_r;

    update_unit #(
        .XY_W_P   (XY_W_P),
        .XY_FRAC_P(XY_FRAC_P),
        .A_BITS_P (A_BITS_P),
        .ACCUM_W_P(ACCUM_W_P),
        .PROD_W_P (PROD_W_P)
    ) u_ph3 (
        .clk     (clk),
        .start   (ph3_start),
        .x_i     (x_i_r),
        .y_i     (y_i_r),
        .Jx_i    (Jx_i_w),
        .a_m     (a_m),
        .dt_fp   (dt_fp),
        .c0_fp   (c0_fp),
        .x_i_new (x_new_w),
        .y_i_new (y_new_w),
        .done    (ph3_done)
    );

    // =========================================================================
    // FSM  (3 states: IDLE → SLC_RUN → UPDATE)
    // =========================================================================
    typedef enum logic [1:0] {
        SC_IDLE    = 2'd0,
        SC_SLC_RUN = 2'd1,
        SC_UPDATE  = 2'd2
    } sc_state_t;

    sc_state_t state;

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state     <= SC_IDLE;
            slc_start <= 1'b0;
            ph3_start <= 1'b0;
            done      <= 1'b0;
            x_i_new   <= '0;
            y_i_new   <= '0;
            x_i_r     <= '0;
            y_i_r     <= '0;
        end else begin
            // default pulse signals
            slc_start <= 1'b0;
            ph3_start <= 1'b0;
            done      <= 1'b0;

            unique case (state)

                // ── Wait for start ────────────────────────────────────────────
                SC_IDLE: begin
                    if (start) begin
                        x_i_r     <= x_i;   // latch inputs
                        y_i_r     <= y_i;
                        slc_start <= 1'b1;  // kick the SLC
                        state     <= SC_SLC_RUN;
                    end
                end

                // ── Wait for SLC + dotprod to finish ──────────────────────────
                // slc_done fires when the last chunk's dotprod is complete
                // and Jx_i_w holds the final result.
                SC_SLC_RUN: begin
                    if (slc_done) begin
                        ph3_start <= 1'b1;  // fire update_unit this cycle
                        state     <= SC_UPDATE;
                    end
                end

                // ── Wait one cycle for update_unit done ───────────────────────
                SC_UPDATE: begin
                    if (ph3_done) begin
                        x_i_new <= x_new_w;
                        y_i_new <= y_new_w;
                        done    <= 1'b1;
                        state   <= SC_IDLE;
                    end
                end

                default: state <= SC_IDLE;

            endcase
        end
    end

endmodule

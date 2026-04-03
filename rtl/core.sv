// =============================================================================
// 7.  dSB Core — 4-cycle pipeline for ONE oscillator i
//
//     Cycle 0 (PRECHARGE) : precharge RBLs; latch sign(x_i)
//     Cycle 1 (PHASE1)    : XNOR fires across full row; capture RBLs
//     Cycle 2 (PHASE2)    : dot-product reconstructed → Jx_i
//     Cycle 3 (PHASE3)    : momentum + position updated; wall applied
// =============================================================================
module core
  import dsb_pkg::*;
#(
  parameter int unsigned N_P       = N,
  parameter int unsigned IC_BITS_P = IC_BITS,
  parameter int unsigned XY_W_P    = XY_W,
  parameter int unsigned XY_FRAC_P = XY_FRAC,
  parameter int unsigned A_BITS_P  = A_BITS,
  parameter int unsigned ACCUM_W_P = ACCUM_W,
  parameter int unsigned PROD_W_P  = PROD_W
)(
  input  logic                        clk,
  input  logic                        rst_n,
  input  logic                        start,
  // J row + neighbour sign write port
  input  logic                        wwl,
  input  logic [N_P*IC_BITS_P-1:0]   wbl_J,
  input  logic [N_P-1:0]              wbl_signs,
  // Current oscillator state
  input  logic signed [XY_W_P-1:0]   x_i,
  input  logic signed [XY_W_P-1:0]   y_i,
  // Schedule and constants
  input  logic [A_BITS_P-1:0]         a_m,
  input  logic [XY_FRAC_P-1:0]        dt_fp,
  input  logic [XY_FRAC_P-1:0]        c0_fp,
  // Outputs
  output logic signed [XY_W_P-1:0]   x_i_new,
  output logic signed [XY_W_P-1:0]   y_i_new,
  output logic                         sign_xi_new,  // sign(x_i_new)
  output logic                         done
);
 
  // sign convention: 1 = positive (spin +1), 0 = negative (spin −1)
  logic sign_xi_cur;
  assign sign_xi_cur  = ~x_i[XY_W_P-1];
  assign sign_xi_new  = ~x_i_new[XY_W_P-1];
 
  // ── Tile interface ─────────────────────────────────────────────────────────
  dsb_tile_if #(.N_P(N_P), .IC_BITS_P(IC_BITS_P)) tile_bus ();
 
  logic precharge_r;
 
  compute_tile #(.N_P(N_P), .IC_BITS_P(IC_BITS_P)) u_tile (
    .clk       (clk),
    .precharge (precharge_r),
    .wwl       (wwl),
    .wbl_J     (wbl_J),
    .wbl_signs (wbl_signs),
    .sign_xi   (sign_xi_cur),
    .tile_bus  (tile_bus.tile_out)
  );
 
  // ── Phase 1 ────────────────────────────────────────────────────────────────
  logic capture_r;
  logic [N_P*IC_BITS_P-1:0] xnor_J;
  logic [N_P-1:0]             sign_eq;
 
  xnor_phase1 #(.N_P(N_P), .IC_BITS_P(IC_BITS_P)) u_ph1 (
    .clk      (clk),
    .capture  (capture_r),
    .tile_bus (tile_bus.ph1_in),
    .xnor_J   (xnor_J),
    .sign_eq  (sign_eq)
  );
 
  // ── Phase 2 ────────────────────────────────────────────────────────────────
  logic                       ph2_start;
  logic signed [ACCUM_W_P-1:0] Jx_i;
  logic                        ph2_done;
  logic                        sign_xi_r;  // latched at start
 
  dotprod_phase2 #(
    .N_P(N_P), .IC_BITS_P(IC_BITS_P), .ACCUM_W_P(ACCUM_W_P)
  ) u_ph2 (
    .clk     (clk),
    .start   (ph2_start),
    .sign_xi (sign_xi_r),
    .xnor_J  (xnor_J),
    .sign_eq (sign_eq),
    .Jx_i    (Jx_i),
    .done    (ph2_done)
  );
 
  // ── Phase 3 ────────────────────────────────────────────────────────────────
  logic                      ph3_start;
  logic signed [XY_W_P-1:0] x_new_w, y_new_w;
  logic                      ph3_done;
  logic signed [XY_W_P-1:0] x_i_r, y_i_r;  // latched at start
 
  update_unit #(
    .XY_W_P(XY_W_P), .XY_FRAC_P(XY_FRAC_P), .A_BITS_P(A_BITS_P),
    .ACCUM_W_P(ACCUM_W_P), .PROD_W_P(PROD_W_P)
  ) u_ph3 (
    .clk    (clk),
    .start  (ph3_start),
    .x_i    (x_i_r),
    .y_i    (y_i_r),
    .Jx_i   (Jx_i),
    .a_m    (a_m),
    .dt_fp  (dt_fp),
    .c0_fp  (c0_fp),
    .x_i_new(x_new_w),
    .y_i_new(y_new_w),
    .done   (ph3_done)
  );
 
  // ── FSM ────────────────────────────────────────────────────────────────────
  core_state_t state;
 
  // ── Combinational control outputs (Moore, from current state) ─────────────
  // Generating these combinatorially rather than as registered NBAs eliminates
  // the silent one-cycle pipeline gap that would otherwise push 'done' one
  // extra cycle beyond the 4-cycle promise (PRECHARGE→PHASE1→PHASE2→PHASE3).
  //
  //   precharge_r : high while FSM is in S_PRECHARGE
  //   capture_r   : high while FSM is in S_PHASE1  (latches RBL on posedge)
  //   ph2_start   : high while FSM is in S_PHASE2  (dotprod FF fires this edge)
  //   ph3_start   : high while FSM is in S_PHASE3 AND ph2_done is already 1
  //                 (update_unit FF fires this edge; FSM will see state→S_IDLE
  //                  at the same posedge via the always_ff below)
  always_comb begin
    precharge_r = (state == S_PRECHARGE);
    capture_r   = (state == S_PHASE1);
    ph2_start   = (state == S_PHASE2);
    ph3_start   = (state == S_PHASE3) && ph2_done;
  end
 
  always_ff @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
      state     <= S_IDLE;
      sign_xi_r <= 1'b0;
      x_i_r     <= '0;
      y_i_r     <= '0;
      x_i_new   <= '0;
      y_i_new   <= '0;
      done      <= 1'b0;
    end else begin
      done <= 1'b0;   // default; overridden by ph3_done branch below
 
      unique case (state)
 
        S_IDLE: begin
          if (start) begin
            x_i_r     <= x_i;
            y_i_r     <= y_i;
            sign_xi_r <= sign_xi_cur;
            state     <= S_PRECHARGE;
          end
        end
 
        S_PRECHARGE: state <= S_PHASE1;
 
        S_PHASE1: state <= S_PHASE2;   // capture_r already high this cycle (comb)
 
        S_PHASE2: state <= S_PHASE3;   // ph2_start already high this cycle (comb)
 
        S_PHASE3: begin
          // ph3_start is combinatorially high this cycle when ph2_done=1;
          // advance to IDLE so the next cycle can catch ph3_done.
          if (ph2_done) state <= S_IDLE;
          // else remain in S_PHASE3 (ph2_done not yet asserted)
        end
 
        default: state <= S_IDLE;
 
      endcase
 
      // Latch phase-3 result whenever it finishes
      if (ph3_done) begin
        x_i_new <= x_new_w;
        y_i_new <= y_new_w;
        done    <= 1'b1;
      end
    end
  end
 
endmodule
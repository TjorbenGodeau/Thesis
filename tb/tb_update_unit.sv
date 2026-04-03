`timescale 1ns / 1ps
// =============================================================================
// TB 6 — dsb_update_unit
//
// Root cause of original failures — timing:
//
//   The always_ff in dsb_update_unit registers x_clamped/y_clamped into
//   x_i_new/y_i_new on the posedge clk where start=1.  Outputs and done
//   are therefore valid AFTER that same rising edge, readable any time
//   before the next rising edge.
//
//   The broken task did:
//     @(posedge clk); #1;    ← past the edge: FF already evaluated start=0
//     start = 1;              ← too late, FF already clocked this cycle
//     @(posedge clk); #1;
//     start = 0;
//     @(posedge clk); #1;   ← reading outputs one cycle after start dropped,
//                              but start was never seen by the FF
//
//   Fix: drive ALL inputs (including start=1) BEFORE the rising edge,
//   then read outputs AFTER the edge where start=1 was sampled.
//   Pattern:
//     1. Set inputs combinatorially (no clock wait needed - comb outputs live)
//     2. Wait for negedge (midpoint) to confirm comb values stable
//     3. @(posedge clk)  ← FF samples start=1 HERE, done/outputs valid after
//     4. Read x_i_new, y_i_new, done  ← all valid now, before next edge
//     5. Drop start=0 for next cycle
//
//   Wall tests (T6.4/T6.5): because x_clamped is purely combinational,
//   we can also verify it directly before the clock edge as a sanity check.
//
// Tests:
//   T6.1  Zero force: Jx=0, a=0, y≠0 → x drifts, y unchanged
//   T6.2  Positive Jx: y grows, x grows
//   T6.3  Negative x with a>0: restoring force adds to Jx drive
//   T6.4  Wall hit positive: x_pre > 1 → x_new=+1.0, y_new=0
//   T6.5  Wall hit negative: x_pre < −1 → x_new=−1.0, y_new=0
//   T6.6  No spurious done without start; done is exactly 1 cycle wide
//   T6.7  Maximum damping: large x with a=a0, Jx=0 → x decays toward 0
// =============================================================================
module tb_update_unit;
  import dsb_pkg::*;
  import tb_helpers::*;
 
  localparam int unsigned XY_T   = XY_W;
  localparam int unsigned FRAC_T = XY_FRAC;
  localparam int unsigned A_T    = A_BITS;
  localparam int unsigned ACC_T  = ACCUM_W;
  localparam int unsigned PROD_T = PROD_W;
 
  localparam int unsigned DT_FP = 1638;   // 0.1 × 2^14
  localparam int unsigned C0_FP = 8192;   // 0.5 × 2^14
  localparam int unsigned A0_FP = 32768;  // 1.0 in Q1.15
 
  logic clk = 0;
  always #5 clk = ~clk;
 
  logic                         start;
  logic signed [XY_T-1:0]      x_i, y_i;
  logic signed [ACC_T-1:0]     Jx_i;
  logic        [A_T-1:0]       a_m;
  logic        [FRAC_T-1:0]    dt_fp, c0_fp;
  logic signed [XY_T-1:0]      x_i_new, y_i_new;
  logic                         done;
 
  update_unit #(
    .XY_W_P(XY_T), .XY_FRAC_P(FRAC_T), .A_BITS_P(A_T),
    .ACCUM_W_P(ACC_T), .PROD_W_P(PROD_T)
  ) dut (.*);
 
  int pass_cnt = 0, fail_cnt = 0;
 
  // ── Reference model ─────────────────────────────────────────────────────
  // Exact mirror of the Python dSB step in real arithmetic.
  task automatic reference_update(
    input  real x_r, y_r, Jx_r, a_r, dt_r, c0_r,
    output real x_exp, y_exp
  );
    real force_val, y_n, x_n;
    force_val = -a_r * x_r + c0_r * Jx_r;
    y_n   = y_r + dt_r * force_val;
    x_n   = x_r + dt_r * y_n;
    if      (x_n >  1.0) begin x_exp =  1.0; y_exp = 0.0; end
    else if (x_n < -1.0) begin x_exp = -1.0; y_exp = 0.0; end
    else                 begin x_exp =  x_n;  y_exp = y_n;  end
  endtask
 
  // ── Drive inputs and pulse start for exactly one clock cycle ─────────────
  // Inputs are set combinatorially (no #delay needed for comb signals).
  // start=1 is in place BEFORE the next posedge so the FF sees it.
  // Outputs are read immediately after that posedge (still before the next).
  task automatic fire_and_check(
    input string label,
    input real   x_r, y_r, Jx_r, a_r
  );
    real x_exp_r, y_exp_r;
    real dt_r, c0_r;
    logic signed [XY_T-1:0] x_exp_fp, y_exp_fp;
 
    dt_r = real'(DT_FP) / real'(1 << FRAC_T);
    c0_r = real'(C0_FP) / real'(1 << FRAC_T);
 
    reference_update(x_r, y_r, Jx_r, a_r, dt_r, c0_r, x_exp_r, y_exp_r);
    x_exp_fp = to_fp(x_exp_r);
    y_exp_fp = to_fp(y_exp_r);
 
    // ── Step 1: set all inputs, assert start ─────────────────────────────
    // Done between clock edges so the FF will sample them on the next posedge.
    @(negedge clk);               // midpoint — safe to drive
    x_i   = to_fp(x_r);
    y_i   = to_fp(y_r);
    Jx_i  = ACC_T'($rtoi(Jx_r));
    a_m   = A_T'($rtoi(a_r * real'(1 << (A_T-1))));
    dt_fp = FRAC_T'(DT_FP);
    c0_fp = FRAC_T'(C0_FP);
    start = 1'b1;
 
    // ── Step 2: rising edge — FF captures inputs, registers update ────────
    @(posedge clk);
    // x_i_new, y_i_new, done are now valid (registered on this edge).
    // Read them after a tiny #1 to be past the NBA region.
    #1;
 
    // ── Step 3: check outputs ─────────────────────────────────────────────
    // Allow ±3 LSB tolerance for fixed-point rounding across multiply chain.
    if ($signed(x_i_new) >= $signed(x_exp_fp) - 3 &&
        $signed(x_i_new) <= $signed(x_exp_fp) + 3 &&
        $signed(y_i_new) >= $signed(y_exp_fp) - 3 &&
        $signed(y_i_new) <= $signed(y_exp_fp) + 3) begin
      $display("  PASS [%s] x=%.4f(exp %.4f)  y=%.4f(exp %.4f)",
               label,
               from_fp(x_i_new), x_exp_r,
               from_fp(y_i_new), y_exp_r);
      pass_cnt++;
    end else begin
      $error("  FAIL [%s] x got=%.4f exp=%.4f  y got=%.4f exp=%.4f",
             label,
             from_fp(x_i_new), x_exp_r,
             from_fp(y_i_new), y_exp_r);
      fail_cnt++;
    end
 
    // ── Step 4: drop start before next posedge ────────────────────────────
    @(negedge clk);
    start = 1'b0;
  endtask
 
  initial begin
    $dumpfile("tb_update_unit.vcd");
    $dumpvars(0, tb_update_unit);
    $display("=== TB6: dsb_update_unit ===");
 
    // Idle defaults
    start = 0; x_i = '0; y_i = '0; Jx_i = '0; a_m = '0;
    dt_fp = FRAC_T'(DT_FP); c0_fp = FRAC_T'(C0_FP);
    repeat(2) @(posedge clk);
 
    // T6.1 Zero force: Jx=0, a=0; y=0.1 so x drifts by dt*y = 0.1*0.1 = 0.01
    fire_and_check("T6.1 zero_force",   0.5,  0.1,  0.0, 0.0);
 
    // T6.2 Positive Jx: force = c0*Jx = 0.5*4 = 2; Δy = dt*2 = 0.2;
    //      y_new=0.2; x_new = 0.1 + dt*0.2 = 0.12
    fire_and_check("T6.2 pos_Jx",       0.1,  0.0,  4.0, 0.0);
 
    // T6.3 Negative x, a=0.5: force = -0.5*(-0.3) + 0.5*2 = 0.15+1 = 1.15
    //      Δy = 0.1*1.15 = 0.115; y_new=0.115; x_new = -0.3+0.0115 = -0.2885
    fire_and_check("T6.3 neg_x",       -0.3,  0.0,  2.0, 0.5);
 
    // T6.4 Wall hit positive: x=0.9, y=5, Jx=4, a=0
    //      force=0.5*4=2; Δy=0.5; y_new=5.5; x_new=0.9+0.55=1.45 → clamp
    //      expected: x=+1.0, y=0
    fire_and_check("T6.4 wall_pos",     0.9,  5.0,  4.0, 0.0);
 
    // T6.5 Wall hit negative: x=-0.9, y=-5, Jx=-4, a=0
    //      force=0.5*(-4)=-2; Δy=-0.2; y_new=-5.2; x_new=-0.9-0.52=-1.42 → clamp
    //      expected: x=-1.0, y=0
    fire_and_check("T6.5 wall_neg",    -0.9, -5.0, -4.0, 0.0);
 
    // ── T6.6a No spurious done ────────────────────────────────────────────
    // start has been 0 for several cycles; done must be 0.
    repeat(2) @(posedge clk); #1;
    if (done === 1'b0) begin
      $display("  PASS [T6.6a no_spurious_done]"); pass_cnt++;
    end else begin
      $error("  FAIL [T6.6a] done=%0b without start", done); fail_cnt++;
    end
 
    // ── T6.6b done is exactly one cycle wide ─────────────────────────────
    // Drive start=1 for one cycle using the same negedge→posedge protocol.
    @(negedge clk);
    x_i = to_fp(0.1); y_i = '0; Jx_i = ACC_T'(4);
    a_m = A_T'(32768); start = 1'b1;
 
    @(posedge clk); #1;          // FF fires: done should be 1 NOW
    if (done === 1'b1) begin
      $display("  PASS [T6.6b done_high_on_start_cycle]"); pass_cnt++;
    end else begin
      $error("  FAIL [T6.6b] done=%0b, expected 1", done); fail_cnt++;
    end
 
    @(negedge clk); start = 1'b0;  // drop start
 
    @(posedge clk); #1;          // next cycle: done must be 0
    if (done === 1'b0) begin
      $display("  PASS [T6.6c done_drops_next_cycle]"); pass_cnt++;
    end else begin
      $error("  FAIL [T6.6c] done still high one cycle after start dropped");
      fail_cnt++;
    end
 
    // T6.7 Maximum damping: x=0.8, y=0, Jx=0, a=1.0
    //      force = -1.0*0.8 = -0.8; Δy = 0.1*(-0.8) = -0.08
    //      y_new = -0.08; x_new = 0.8 + 0.1*(-0.08) = 0.792
    fire_and_check("T6.7 max_damping",  0.8,  0.0,  0.0, 1.0);
 
    $display("--- TB6 result: %0d PASS, %0d FAIL ---\n", pass_cnt, fail_cnt);
    #20 $finish;
  end
 
endmodule : tb_update_unit
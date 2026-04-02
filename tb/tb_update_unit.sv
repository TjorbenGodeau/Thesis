`timescale 1ns / 1ps
// =============================================================================
// TB 6 — dsb_update_unit
// Tests verify the full dSB momentum/position update in fixed-point.
// Reference model computes the same result in real arithmetic and converts.
//
//   T6.1  Zero force: Jx=0, a=0 → x, y unchanged
//   T6.2  Positive Jx only: y grows, x grows
//   T6.3  Negative x with positive Jx: force = c0*Jx - a*x (a*(-x) positive)
//   T6.4  Wall hit positive: x_pre > 1 → x_new=+1, y_new=0
//   T6.5  Wall hit negative: x_pre < −1 → x_new=−1, y_new=0
//   T6.6  done pulses correctly
//   T6.7  a=a0 (maximum): damping term dominates for large x
// =============================================================================
module tb_update_unit;
  import dsb_pkg::*;
  import tb_helpers::*;
 
  localparam int unsigned XY_T   = XY_W;
  localparam int unsigned FRAC_T = XY_FRAC;
  localparam int unsigned A_T    = A_BITS;
  localparam int unsigned ACC_T  = ACCUM_W;
  localparam int unsigned PROD_T = PROD_W;
 
  // Fixed-point constants matching default dSB
  localparam int unsigned DT_FP = 1638;   // 0.1  × 2^14
  localparam int unsigned C0_FP = 8192;   // 0.5  × 2^14
  localparam int unsigned A0_FP = 32768;  // 1.0  in Q1.15
 
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
    .XY_W_P(XY_T), 
    .XY_FRAC_P(FRAC_T), 
    .A_BITS_P(A_T),
    .ACCUM_W_P(ACC_T), 
    .PROD_W_P(PROD_T)
  ) dut (.*);
 
  int pass_cnt = 0, fail_cnt = 0;
 
  // ── Reference model in real arithmetic ────────────────────────────────────
  // Mirrors exactly the Python dSB step:
  //   a_real = a_m / 2^(A_BITS-1)
  //   force  = -a_real * x_real + c0_real * Jx_real
  //   y_new  = y_real + dt_real * force
  //   x_new  = x_real + dt_real * y_new
  //   wall clamp
  task automatic reference_update(
    input  real x_r, y_r, Jx_r,
    input  real a_r, dt_r, c0_r,
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
 
    @(posedge clk); #1;
    x_i   = to_fp(x_r);
    y_i   = to_fp(y_r);
    Jx_i  = ACC_T'($rtoi(Jx_r));
    a_m   = A_T'($rtoi(a_r * real'(1 << (A_T-1))));
    dt_fp = FRAC_T'(DT_FP);
    c0_fp = FRAC_T'(C0_FP);
    start = 1;
    @(posedge clk); #1;
    start = 0;
    @(posedge clk); #1;   // wait 1 cycle for done
 
    // Allow ±2 LSB tolerance for fixed-point rounding
    if ($signed(x_i_new) >= $signed(x_exp_fp) - 2 &&
        $signed(x_i_new) <= $signed(x_exp_fp) + 2 &&
        $signed(y_i_new) >= $signed(y_exp_fp) - 2 &&
        $signed(y_i_new) <= $signed(y_exp_fp) + 2) begin
      $display("  PASS [%s] x=%.4f(%.4f) y=%.4f(%.4f)",
        label, from_fp(x_i_new), x_exp_r, from_fp(y_i_new), y_exp_r);
      pass_cnt++;
    end else begin
      $error("  FAIL [%s] x got=%.4f exp=%.4f  y got=%.4f exp=%.4f",
        label, from_fp(x_i_new), x_exp_r, from_fp(y_i_new), y_exp_r);
      fail_cnt++;
    end
  endtask
 
  initial begin
    $dumpfile("tb_update_unit.vcd");
    $dumpvars(0, tb_update_unit);
    $display("=== TB6: dsb_update_unit ===");
 
    start = 0; x_i = 0; y_i = 0; Jx_i = 0; a_m = 0;
    dt_fp = FRAC_T'(DT_FP); c0_fp = FRAC_T'(C0_FP);
    repeat(2) @(posedge clk);
 
    // T6.1 Zero force: x=0.5, y=0.1, Jx=0, a=0 → force=0, y unchanged, x grows
    fire_and_check("T6.1 zero_force",  0.5,  0.1, 0.0, 0.0);
 
    // T6.2 Positive Jx drives y positive
    fire_and_check("T6.2 pos_Jx",      0.1,  0.0, 4.0, 0.0);
 
    // T6.3 Negative x: damping term −a*x is positive when x<0
    fire_and_check("T6.3 neg_x",      -0.3,  0.0, 2.0, 0.5);
 
    // T6.4 Wall hit positive: large y forces x past +1
    fire_and_check("T6.4 wall_pos",    0.9,  5.0, 4.0, 0.0);
 
    // T6.5 Wall hit negative
    fire_and_check("T6.5 wall_neg",   -0.9, -5.0, -4.0, 0.0);
 
    // T6.6 done pulse check
    start = 0; repeat(2) @(posedge clk);
    if (done === 1'b0) begin
      $display("  PASS [T6.6 no_spurious_done]"); pass_cnt++;
    end else begin
      $error("  FAIL [T6.6] done high without start"); fail_cnt++;
    end
    @(posedge clk); #1; x_i = to_fp(0.1); y_i = '0;
    Jx_i = 4; a_m = A_T'(32768); start = 1;
    @(posedge clk); #1; start = 0;
    @(posedge clk); #1;
    if (done === 1'b1) begin
      $display("  PASS [T6.6 done_pulses]"); pass_cnt++;
    end else begin
      $error("  FAIL [T6.6 done_pulses] done=%0b", done); fail_cnt++;
    end
 
    // T6.7 Maximum a (full damping), large x → x should decrease
    fire_and_check("T6.7 max_damping",  0.8,  0.0, 0.0, 1.0);
 
    $display("--- TB6 result: %0d PASS, %0d FAIL ---\n", pass_cnt, fail_cnt);
    #20 $finish;
  end
 
endmodule : tb_update_unit
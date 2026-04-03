`timescale 1ns / 1ps
// =============================================================================
// TB 7 — dsb_core
// Tests the 4-cycle FSM for a single oscillator end-to-end.
//
//   T7.1  Single update: load tile, run → verify x, y change correctly
//   T7.2  sign_xi_new tracks sign of x_i_new
//   T7.3  Wall: start near boundary → x clamps to ±1
//   T7.4  Sequential updates: run two consecutive updates, state evolves
//   T7.5  done is only high for exactly 1 cycle after phase 3
// =============================================================================
module tb_core;
  import dsb_pkg::*;
  import tb_helpers::*;
 
  localparam int unsigned N_T   = 4;
  localparam int unsigned IC_T  = 8;
  localparam int unsigned XY_T  = XY_W;
  localparam int unsigned FT    = XY_FRAC;
  localparam int unsigned AT    = A_BITS;
  localparam int unsigned ACCT  = ACCUM_W;
  localparam int unsigned PRODT = PROD_W;
 
  localparam int unsigned DT_FP  = 1638;
  localparam int unsigned C0_FP  = 8192;
  localparam int unsigned A0_FP  = 32768;
 
  logic clk = 0, rst_n;
  always #5 clk = ~clk;
 
  logic                       start, wwl;
  logic [N_T*IC_T-1:0]       wbl_J;
  logic [N_T-1:0]             wbl_signs;
  logic signed [XY_T-1:0]    x_i, y_i;
  logic [AT-1:0]              a_m;
  logic [FT-1:0]              dt_fp, c0_fp;
  logic signed [XY_T-1:0]    x_i_new, y_i_new;
  logic                       sign_xi_new, done;
 
  core #(
    .N_P(N_T), .IC_BITS_P(IC_T), .XY_W_P(XY_T), .XY_FRAC_P(FT),
    .A_BITS_P(AT), .ACCUM_W_P(ACCT), .PROD_W_P(PRODT)
  ) dut (.*);
 
  int pass_cnt = 0, fail_cnt = 0;
 
  task automatic reset_dut;
    rst_n = 0; start = 0; wwl = 0;
    wbl_J = '0; wbl_signs = '0; x_i = '0; y_i = '0;
    a_m = '0; dt_fp = FT'(DT_FP); c0_fp = FT'(C0_FP);
    repeat(3) @(posedge clk);
    rst_n = 1;
    repeat(2) @(posedge clk);
  endtask
 
  // Load J row and sign snapshot into tile, then run one update
  task automatic run_one(
    input logic [N_T*IC_T-1:0] j_row,
    input logic [N_T-1:0]       signs,
    input logic signed [XY_T-1:0] x_val, y_val,
    input logic [AT-1:0]         a_val
  );
    // Write tile
    @(posedge clk); #1;
    wwl = 1; wbl_J = j_row; wbl_signs = signs;
    @(posedge clk); #1; wwl = 0;
    // Set oscillator state and schedule
    x_i = x_val; y_i = y_val; a_m = a_val;
    // Start core
    @(posedge clk); #1; start = 1;
    @(posedge clk); #1; start = 0;
    // Wait for done (max 6 cycles)
    repeat(6) begin
      @(posedge clk); #1;
      if (done) return;
    end
    $error("  FAIL: core did not assert done within 6 cycles");
    fail_cnt++;
  endtask
 
  initial begin
    $dumpfile("tb_core.vcd");
    $dumpvars(0, tb_core);
    $display("=== TB7: dsb_core ===");
    reset_dut();
 
    // T7.1 Single update: x=0.1, y=0, all J=+1, all signs=1 (same as xi),
    //      a=0 → Jx=4×1=4; force = c0×4 = 2; Δy = 0.1×2 = 0.2; y_new=0.2
    //      x_new = 0.1 + 0.1×0.2 = 0.12
    run_one({4{8'h01}}, 4'b1111, to_fp(0.1), '0, '0);
    begin
      real x_got;
      real y_got;
      x_got = from_fp(x_i_new);
      y_got = from_fp(y_i_new);
      if (x_got > 0.10 && x_got < 0.14 && y_got > 0.15 && y_got < 0.25) begin
        $display("  PASS [T7.1] x=%.4f y=%.4f", x_got, y_got); pass_cnt++;
      end else begin
        $error("  FAIL [T7.1] x=%.4f (exp~0.12) y=%.4f (exp~0.20)", x_got, y_got);
        fail_cnt++;
      end
    end
 
    // T7.2 sign_xi_new: positive x_new → sign=1
    if (sign_xi_new === 1'b1) begin
      $display("  PASS [T7.2 sign_pos]"); pass_cnt++;
    end else begin
      $error("  FAIL [T7.2 sign_pos] sign_xi_new=%0b", sign_xi_new); fail_cnt++;
    end
 
    // T7.3 Wall hit: x near +1, large positive Jx, a=0
    //      Jx=4, c0=0.5 → force=2.0; Δy=dt×force=0.1×2.0=0.2; y_new=0.2
    //      x_new = 0.99 + 0.1×0.2 = 1.01 > 1.0 → wall clamp → x=+1.0, y=0
    run_one({4{8'h01}}, 4'b1111, to_fp(0.99), '0, '0);
    begin
      real x_got;
      real y_got;
      x_got = from_fp(x_i_new);
      y_got = from_fp(y_i_new);
      if (x_got >= 0.999 && y_got == 0.0) begin
        $display("  PASS [T7.3 wall_pos] x=%.4f y=%.4f", x_got, y_got); pass_cnt++;
      end else begin
        $error("  FAIL [T7.3 wall_pos] x=%.4f y=%.4f", x_got, y_got); fail_cnt++;
      end
    end
 
    // T7.4 Sequential updates: run twice from same starting point and verify
    //      second update uses first output
    reset_dut();
    run_one({4{8'h01}}, 4'b1111, to_fp(0.05), '0, '0);
    begin
      logic signed [XY_T-1:0] x1;
      logic signed [XY_T-1:0] y1;
      real x2;

      x1 = x_i_new;
      y1 = y_i_new;
      run_one({4{8'h01}}, 4'b1111, x1, y1, '0);
      x2 = from_fp(x_i_new);
      if (x2 > from_fp(x1)) begin
        $display("  PASS [T7.4 sequential] x1=%.4f x2=%.4f", from_fp(x1), x2);
        pass_cnt++;
      end else begin
        $error("  FAIL [T7.4 sequential] x should grow but x1=%.4f x2=%.4f",
               from_fp(x1), x2);
        fail_cnt++;
      end
    end
 
    // T7.5 done is a single-cycle pulse
    reset_dut();
    run_one({4{8'h01}}, 4'b1111, to_fp(0.1), '0, '0);
    // done was high last cycle; check it drops the next cycle
    @(posedge clk); #1;
    if (done === 1'b0) begin
      $display("  PASS [T7.5 done_single_cycle]"); pass_cnt++;
    end else begin
      $error("  FAIL [T7.5 done_single_cycle] done still high"); fail_cnt++;
    end
 
    $display("--- TB7 result: %0d PASS, %0d FAIL ---\n", pass_cnt, fail_cnt);
    #20 $finish;
  end
 
endmodule : tb_core
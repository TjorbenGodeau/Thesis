`timescale 1ns / 1ps
// =============================================================================
// TB 3 — dsb_xnor_phase1
// Tests:
//   T3.1  capture=0 → outputs hold previous value (no latch)
//   T3.2  capture=1 → outputs updated to tile_bus values
//   T3.3  Burst: two different tile_bus values, only second captured
// =============================================================================
module tb_xnor_phase1;
  import dsb_pkg::*;
 
  localparam int unsigned N_T  = 4;
  localparam int unsigned IC_T = 4;
 
  logic clk = 0;
  always #5 clk = ~clk;
 
  logic capture;
  dsb_tile_if #(.N_P(N_T), .IC_BITS_P(IC_T)) tile_bus ();
 
  logic [N_T*IC_T-1:0] xnor_J;
  logic [N_T-1:0]       sign_eq;
 
  xnor_phase1 #(.N_P(N_T), .IC_BITS_P(IC_T)) dut (
    .clk      (clk),
    .capture  (capture),
    .tile_bus (tile_bus.ph1_in),
    .xnor_J   (xnor_J),
    .sign_eq  (sign_eq)
  );

  logic [N_T*IC_T-1:0] tb_rbl_J    = '0;
  logic [N_T-1:0]       tb_rbl_signs = '0;
  
  // Drive tile_bus directly from TB (bypass tile module)
  assign tile_bus.rbl_J    = tb_rbl_J;
  assign tile_bus.rbl_signs = tb_rbl_signs;
 
  int pass_cnt = 0, fail_cnt = 0;
 
  task automatic check_capture(
    input string             label,
    input logic [N_T*IC_T-1:0] got_J,  exp_J,
    input logic [N_T-1:0]      got_s,  exp_s
  );
    if (got_J === exp_J && got_s === exp_s) begin
      $display("  PASS [%s]", label); pass_cnt++;
    end else begin
      $error("  FAIL [%s] xnor_J got=%0b exp=%0b  sign_eq got=%0b exp=%0b",
             label, got_J, exp_J, got_s, exp_s);
      fail_cnt++;
    end
  endtask
 
  initial begin
    $dumpfile("tb_xnor_phase1.vcd");
    $dumpvars(0, tb_xnor_phase1);
    $display("=== TB3: dsb_xnor_phase1 ===");
 
    capture = 0;
    repeat(2) @(posedge clk);
 
    // T3.1 capture=0: drive new value, output must not change
    tb_rbl_J = {N_T*IC_T{1'b1}}; tb_rbl_signs = {N_T{1'b1}};
    @(posedge clk); #1;
    check_capture("T3.1 no_capture",
                  xnor_J, {N_T*IC_T{1'b0}},
                  sign_eq, {N_T{1'b0}});
 
    // T3.2 capture=1: output should latch
    @(posedge clk); #1; capture = 1;
    @(posedge clk); #1; capture = 0;
    #1;
    check_capture("T3.2 capture_ones",
                  xnor_J, {N_T*IC_T{1'b1}},
                  sign_eq, {N_T{1'b1}});
 
    // T3.3 Change bus without capture, then capture only second value
    tb_rbl_J = {N_T*IC_T{1'b0}}; tb_rbl_signs = {N_T{1'b0}};
    @(posedge clk); #1;   // no capture
    tb_rbl_J = 16'hA5A5;  tb_rbl_signs = 4'b1100;
    @(posedge clk); #1; capture = 1;
    @(posedge clk); #1; capture = 0;
    #1;
    check_capture("T3.3 capture_second",
                  xnor_J, 16'hA5A5,
                  sign_eq, 4'b1100);
 
    $display("--- TB3 result: %0d PASS, %0d FAIL ---\n", pass_cnt, fail_cnt);
    #20 $finish;
  end
 
endmodule : tb_xnor_phase1
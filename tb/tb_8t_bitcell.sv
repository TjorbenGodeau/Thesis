`timescale 1ns / 1ps
// =============================================================================
// TB 1 — dsb_8t_bitcell
// Tests:
//   T1.1  Write 0, compute XNOR with rwl=0 → expect 1 (match)
//   T1.2  Write 0, compute XNOR with rwl=1 → expect 0 (mismatch)
//   T1.3  Write 1, compute XNOR with rwl=1 → expect 1 (match)
//   T1.4  Write 1, compute XNOR with rwl=0 → expect 0 (mismatch)
//   T1.5  Precharge overrides compute → expect 1 regardless of rwl
//   T1.6  wwl=1 during compute phase → rbl held high (normal mode)
//   T1.7  Write sequence: overwrite 0→1, verify new stored value
// =============================================================================
module tb_8t_bitcell;
  import dsb_pkg::*;
 
  logic clk = 0;
  always #5 clk = ~clk;
 
  logic wwl, wbl, rwl, precharge;
  logic rbl_out;
 
  dsb_8t_bitcell dut (.*);
 
  // ── Pass/fail counter ──────────────────────────────────────────────────────
  int pass_cnt = 0, fail_cnt = 0;
 
  task automatic check(
    input string test_name,
    input logic  got,
    input logic  exp
  );
    if (got === exp) begin
      $display("  PASS [%s] got=%0b exp=%0b", test_name, got, exp);
      pass_cnt++;
    end else begin
      $error("  FAIL [%s] got=%0b exp=%0b", test_name, got, exp);
      fail_cnt++;
    end
  endtask
 
  // ── Write helper: drive wwl+wbl for one cycle ─────────────────────────────
  task automatic write_cell(input logic val);
    @(posedge clk); #1;
    wwl = 1; wbl = val; rwl = 0; precharge = 0;
    @(posedge clk); #1;
    wwl = 0; wbl = 0;
  endtask
 
  initial begin
    $dumpfile("tb_8t_bitcell.vcd");
    $dumpvars(0, tb_8t_bitcell);
    $display("=== TB1: dsb_8t_bitcell ===");
 
    // Default idle state
    wwl = 0; wbl = 0; rwl = 0; precharge = 0;
    repeat(2) @(posedge clk);
 
    // T1.1 Write 0, XNOR rwl=0 → match → 1
    write_cell(0);
    @(posedge clk); #1; precharge = 0; wwl = 0; rwl = 0;
    #1; check("T1.1 W0_R0_match",  rbl_out, 1'b1);
 
    // T1.2 Write 0, XNOR rwl=1 → mismatch → 0
    @(posedge clk); #1; rwl = 1;
    #1; check("T1.2 W0_R1_miss",   rbl_out, 1'b0);
 
    // T1.3 Write 1, XNOR rwl=1 → match → 1
    write_cell(1);
    @(posedge clk); #1; precharge = 0; wwl = 0; rwl = 1;
    #1; check("T1.3 W1_R1_match",  rbl_out, 1'b1);
 
    // T1.4 Write 1, XNOR rwl=0 → mismatch → 0
    @(posedge clk); #1; rwl = 0;
    #1; check("T1.4 W1_R0_miss",   rbl_out, 1'b0);
 
    // T1.5 Precharge forces high regardless of stored value or rwl
    @(posedge clk); #1; precharge = 1; rwl = 0;
    #1; check("T1.5 precharge_rwl0", rbl_out, 1'b1);
    rwl = 1;
    #1; check("T1.5 precharge_rwl1", rbl_out, 1'b1);
    precharge = 0;
 
    // T1.6 wwl=1 holds rbl high during normal mode (no compute discharge)
    @(posedge clk); #1; wwl = 1; wbl = 0; rwl = 1;
    #1; check("T1.6 wwl_high_hold",  rbl_out, 1'b1);
    @(posedge clk); #1; wwl = 0;
 
    // T1.7 Overwrite 0 → 1, verify compute sees new value
    write_cell(0);
    write_cell(1);   // overwrite with 1
    @(posedge clk); #1; wwl = 0; rwl = 1;
    #1; check("T1.7 overwrite_1_R1", rbl_out, 1'b1);
    rwl = 0;
    #1; check("T1.7 overwrite_1_R0", rbl_out, 1'b0);
 
    $display("--- TB1 result: %0d PASS, %0d FAIL ---\n", pass_cnt, fail_cnt);
    #20 $finish;
  end
 
endmodule : tb_8t_bitcell
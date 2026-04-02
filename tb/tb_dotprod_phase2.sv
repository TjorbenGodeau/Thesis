// =============================================================================
// TB 4 — dsb_dotprod_phase2
// Tests (N=4, IC_BITS=8):
//   T4.1  All J=+1 (8'h01), all sign_eq=1, sign_xi=1 → Jx = +4
//   T4.2  All J=+1, all sign_eq=0, sign_xi=1 → all products negated → Jx = -4
//   T4.3  All J=+1, all sign_eq=1, sign_xi=0 → same magnitude, +N correction
//   T4.4  Mixed: J=[+1,−1,+1,−1], sign_eq=1111 → Jx = 0
//   T4.5  done pulses only when start=1
//   T4.6  J=+2 all, sign_eq=1111, sign_xi=1 → Jx = +8
// =============================================================================
module tb_dotprod_phase2;
  import dsb_pkg::*;
 
  localparam int unsigned N_T      = 4;
  localparam int unsigned IC_T     = 8;
  localparam int unsigned ACCUM_T  = 24;
 
  logic clk = 0;
  always #5 clk = ~clk;
 
  logic                      start, sign_xi;
  logic [N_T*IC_T-1:0]      xnor_J;
  logic [N_T-1:0]            sign_eq;
  logic signed [ACCUM_T-1:0] Jx_i;
  logic                       done;
 
  dsb_dotprod_phase2 #(.N_P(N_T), .IC_BITS_P(IC_T), .ACCUM_W_P(ACCUM_T)) dut (.*);
 
  int pass_cnt = 0, fail_cnt = 0;
 
  task automatic fire(
    input logic [N_T*IC_T-1:0] j_in,
    input logic [N_T-1:0]       seq_in,
    input logic                  sxi
  );
    @(posedge clk); #1;
    xnor_J  = j_in;
    sign_eq = seq_in;
    sign_xi = sxi;
    start   = 1;
    @(posedge clk); #1;
    start = 0;
    @(posedge clk); #1;   // wait for done
  endtask
 
  task automatic check_jx(
    input string              label,
    input logic signed [ACCUM_T-1:0] got,
    input int                 exp
  );
    if ($signed(got) === exp) begin
      $display("  PASS [%s] Jx=%0d", label, $signed(got)); pass_cnt++;
    end else begin
      $error("  FAIL [%s] Jx got=%0d exp=%0d", label, $signed(got), exp);
      fail_cnt++;
    end
  endtask
 
  initial begin
    $dumpfile("tb_dotprod_phase2.vcd");
    $dumpvars(0, tb_dotprod_phase2);
    $display("=== TB4: dsb_dotprod_phase2 ===");
 
    start = 0; sign_xi = 0; xnor_J = '0; sign_eq = '0;
    repeat(2) @(posedge clk);
 
    // T4.1 All J=+1 (XNOR result = 8'h01), all same sign, xi positive
    // Each corrected word = 8'h01; sign-extended = +1; sum = 4; no correction
    fire({4{8'h01}}, 4'b1111, 1'b1);
    check_jx("T4.1 all_plus1_seq1_xi1", Jx_i, 4);
 
    // T4.2 All J=+1, sign_eq=0 (diff sign) → corrected = ~8'h01 = 8'hFE = −2 signed
    // sum = 4×(−2) = −8; no correction (xi=1)
    fire({4{8'h01}}, 4'b0000, 1'b1);
    check_jx("T4.2 all_plus1_seq0_xi1", Jx_i, -8);
 
    // T4.3 All J=+1, sign_eq=1111, xi=0 → same as T4.1 but +N added
    // result = 4 + 4 = 8
    fire({4{8'h01}}, 4'b1111, 1'b0);
    check_jx("T4.3 all_plus1_seq1_xi0", Jx_i, 8);
 
    // T4.4 Mixed: J=[+1,−1,+1,−1] as [8'h01,8'hFF,8'h01,8'hFF], seq=1111, xi=1
    // signed values: +1 + (−1) + (+1) + (−1) = 0
    fire({8'hFF, 8'h01, 8'hFF, 8'h01}, 4'b1111, 1'b1);
    check_jx("T4.4 mixed_pm1", Jx_i, 0);
 
    // T4.5 done should only pulse when start fires
    start = 0;
    repeat(3) @(posedge clk);
    if (done === 1'b0) begin
      $display("  PASS [T4.5 no_spurious_done]"); pass_cnt++;
    end else begin
      $error("  FAIL [T4.5 no_spurious_done] done was high without start");
      fail_cnt++;
    end
 
    // T4.6 J=+2 (8'h02), seq=1111, xi=1 → 4×2 = 8
    fire({4{8'h02}}, 4'b1111, 1'b1);
    check_jx("T4.6 all_plus2", Jx_i, 8);
 
    $display("--- TB4 result: %0d PASS, %0d FAIL ---\n", pass_cnt, fail_cnt);
    #20 $finish;
  end
 
endmodule : tb_dotprod_phase2
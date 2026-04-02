`timescale 1ns / 1ps
// =============================================================================
// TB 5 — dsb_schedule
// Tests:
//   T5.1  m=0         → a_m = a0_fp  (maximum)
//   T5.2  m=Nstep/2   → a_m ≈ a0_fp/2
//   T5.3  m=Nstep-1   → a_m ≈ a0_fp/Nstep (one step remaining)
//   T5.4  m=Nstep     → a_m = 0  (schedule complete)
//   T5.5  m > Nstep   → a_m = 0  (clamp)
//   T5.6  Monotone decrease across full schedule
// =============================================================================
module tb_schedule;
  import dsb_pkg::*;
 
  localparam int unsigned A_T    = 16;
  localparam int unsigned STEP_T = 12;
  localparam int unsigned NSTEP  = 50;
  localparam int unsigned A0_FP  = 32768;  // 1.0 in Q1.15
 
  logic [STEP_T-1:0] m;
  logic [STEP_T-1:0] Nstep_in = STEP_T'(NSTEP);
  logic [A_T-1:0]    a0_fp_in = A_T'(A0_FP);
  logic [A_T-1:0]    a_m;
 
  schedule #(.A_BITS_P(A_T), .STEP_W_P(STEP_T)) dut (
    .m     (m),
    .Nstep (Nstep_in),
    .a0_fp (a0_fp_in),
    .a_m   (a_m)
  );
 
  int pass_cnt = 0, fail_cnt = 0;
 
  // Expected value: a0_fp * (Nstep - m) / Nstep (integer division)
  function automatic int unsigned expected_a(input int unsigned m_v);
    return (A0_FP * (NSTEP - m_v)) / NSTEP;
  endfunction
 
  task automatic check_a(input string label, input int unsigned exp);
    if (a_m === A_T'(exp)) begin
      $display("  PASS [%s] a_m=%0d (exp %0d)", label, a_m, exp); pass_cnt++;
    end else begin
      // Allow ±1 LSB for integer division rounding
      if (a_m >= A_T'(exp-1) && a_m <= A_T'(exp+1)) begin
        $display("  PASS [%s] a_m=%0d (exp %0d, within ±1 LSB)", label, a_m, exp);
        pass_cnt++;
      end else begin
        $error("  FAIL [%s] a_m=%0d exp=%0d", label, a_m, exp);
        fail_cnt++;
      end
    end
  endtask
 
  initial begin
    $dumpfile("tb_schedule.vcd");
    $dumpvars(0, tb_schedule);
    $display("=== TB5: dsb_schedule ===");
 
    // T5.1 m=0 → maximum a
    m = 0; #1;
    check_a("T5.1 m=0",      expected_a(0));
 
    // T5.2 m=Nstep/2
    m = STEP_T'(NSTEP/2); #1;
    check_a("T5.2 m=half",   expected_a(NSTEP/2));
 
    // T5.3 m=Nstep-1
    m = STEP_T'(NSTEP-1); #1;
    check_a("T5.3 m=Nstep-1", expected_a(NSTEP-1));
 
    // T5.4 m=Nstep → 0
    m = STEP_T'(NSTEP); #1;
    if (a_m === '0) begin
      $display("  PASS [T5.4 m=Nstep → 0]"); pass_cnt++;
    end else begin
      $error("  FAIL [T5.4 m=Nstep] a_m=%0d exp=0", a_m); fail_cnt++;
    end
 
    // T5.5 m > Nstep → still 0
    m = STEP_T'(NSTEP + 5); #1;
    if (a_m === '0) begin
      $display("  PASS [T5.5 m>Nstep → 0]"); pass_cnt++;
    end else begin
      $error("  FAIL [T5.5 m>Nstep] a_m=%0d exp=0", a_m); fail_cnt++;
    end
 
    // T5.6 Monotone decrease
    begin
      int unsigned prev = A0_FP + 1;
      int unsigned curr;
      logic ok = 1;
      for (int i = 0; i <= NSTEP; i++) begin
        m = STEP_T'(i); #1;
        curr = unsigned'(a_m);
        if (curr > prev) begin ok = 0; break; end
        prev = curr;
      end
      if (ok) begin
        $display("  PASS [T5.6 monotone_decrease]"); pass_cnt++;
      end else begin
        $error("  FAIL [T5.6 monotone_decrease] not monotonically decreasing");
        fail_cnt++;
      end
    end
 
    $display("--- TB5 result: %0d PASS, %0d FAIL ---\n", pass_cnt, fail_cnt);
    #20 $finish;
  end
 
endmodule : tb_schedule
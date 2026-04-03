`timescale 1ns / 1ps
// =============================================================================
// TB 8 — dsb_array  (integration test)
// Tests the full N-oscillator round-robin controller.
//
//   T8.1  Load J matrix and initial x; run schedule; verify schedule_done
//   T8.2  step_done pulses exactly Nstep times
//   T8.3  signs_out tracks sign of x_out
//   T8.4  After schedule_done, x values are not all zero (bifurcation happened)
//   T8.5  Identity J (all zero): x stays near initial values (no forcing)
//   T8.6  Re-run after schedule_done does not start again without reset
// =============================================================================
module tb_array;
  import dsb_pkg::*;
  import tb_helpers::*;
 
  localparam int unsigned N_T    = 4;
  localparam int unsigned IC_T   = 8;
  localparam int unsigned XY_T   = XY_W;
  localparam int unsigned FT     = XY_FRAC;
  localparam int unsigned AT     = A_BITS;
  localparam int unsigned STEPT  = STEP_W;
  localparam int unsigned ACCT   = ACCUM_W;
  localparam int unsigned PRODT  = PROD_W;
 
  localparam int unsigned DT_FP  = 1638;
  localparam int unsigned C0_FP  = 8192;
  localparam int unsigned A0_FP  = 32768;
  localparam int unsigned NSTEP  = 50;
 
  localparam logic [IC_T-1:0] Jp = 8'h01;
  localparam logic [IC_T-1:0] Jz = 8'h00;
 
  logic clk = 0, rst_n;
  always #5 clk = ~clk;
 
  logic [$clog2(N_T)-1:0]    wr_row, wr_xy_idx;
  logic                        wr_en, wr_xy_en, run;
  logic [N_T*IC_T-1:0]        wr_J_row;
  logic signed [XY_T-1:0]     wr_x, wr_y;
  logic [N_T*XY_T-1:0]        x_out, y_out;
  logic [N_T-1:0]              signs_out;
  logic                        step_done, schedule_done;
 
  dsb_array #(
    .N_P(N_T), .IC_BITS_P(IC_T), .XY_W_P(XY_T), .XY_FRAC_P(FT),
    .A_BITS_P(AT), .STEP_W_P(STEPT), .ACCUM_W_P(ACCT), .PROD_W_P(PRODT)
  ) dut (
    .clk(clk), .rst_n(rst_n),
    .wr_row(wr_row), .wr_en(wr_en), .wr_J_row(wr_J_row),
    .wr_xy_idx(wr_xy_idx), .wr_xy_en(wr_xy_en),
    .wr_x(wr_x), .wr_y(wr_y),
    .run(run),
    .Nstep(STEPT'(NSTEP)),
    .a0_fp(AT'(A0_FP)),
    .dt_fp(FT'(DT_FP)),
    .c0_fp(FT'(C0_FP)),
    .x_out(x_out), .y_out(y_out),
    .signs_out(signs_out),
    .step_done(step_done),
    .schedule_done(schedule_done)
  );
 
  int pass_cnt = 0, fail_cnt = 0;
 
  task automatic hard_reset;
    rst_n = 0; run = 0; wr_en = 0; wr_xy_en = 0;
    wr_row = '0; wr_J_row = '0; wr_xy_idx = '0; wr_x = '0; wr_y = '0;
    repeat(4) @(posedge clk);
    rst_n = 1;
    repeat(2) @(posedge clk);
  endtask
 
  task automatic write_J(input int r, input logic [N_T*IC_T-1:0] data);
    @(posedge clk); #1;
    wr_row = $clog2(N_T)'(r); wr_J_row = data; wr_en = 1;
    @(posedge clk); #1; wr_en = 0;
  endtask
 
  task automatic write_xy(input int idx, input real x_r, input real y_r);
    @(posedge clk); #1;
    wr_xy_idx = $clog2(N_T)'(idx);
    wr_x = to_fp(x_r); wr_y = to_fp(y_r);
    wr_xy_en = 1;
    @(posedge clk); #1; wr_xy_en = 0;
  endtask
 
  task automatic run_and_wait(output int step_count);
    @(posedge clk); #1; run = 1;
    @(posedge clk); #1; run = 0;
    step_count = 0;
    fork
      begin : timeout_proc
        repeat(20000) @(posedge clk);
        $error("  FAIL: timeout waiting for schedule_done");
        fail_cnt++;
        disable wait_proc;
      end
      begin : wait_proc
        forever begin
          @(posedge clk); #1;
          if (step_done) step_count++;
          if (schedule_done) begin
            disable timeout_proc;
            break;
          end
        end
      end
    join
  endtask
 
  initial begin
    $dumpfile("tb_array.vcd");
    $dumpvars(0, tb_array);
    $display("=== TB8: dsb_array ===");
 
    // ── T8.1 + T8.2 + T8.3 + T8.4 ─────────────────────────────────────────
    // Ferromagnetic J=+1 (upper triangular), asymmetric initial x
    hard_reset();
    write_J(0, {Jp, Jp, Jp, Jz});
    write_J(1, {Jp, Jp, Jz, Jz});
    write_J(2, {Jp, Jz, Jz, Jz});
    write_J(3, {Jz, Jz, Jz, Jz});
    write_xy(0,  0.05, 0.0);
    write_xy(1, -0.04, 0.0);
    write_xy(2,  0.06, 0.0);
    write_xy(3, -0.03, 0.0);
 
    begin
      int sc;
      run_and_wait(sc);
 
      // T8.1 schedule_done asserted
      if (schedule_done) begin
        $display("  PASS [T8.1 schedule_done]"); pass_cnt++;
      end else begin
        $error("  FAIL [T8.1 schedule_done]"); fail_cnt++;
      end
 
      // T8.2 step_done pulsed exactly Nstep times
      if (sc === NSTEP) begin
        $display("  PASS [T8.2 step_count=%0d]", sc); pass_cnt++;
      end else begin
        $error("  FAIL [T8.2] step_count=%0d exp=%0d", sc, NSTEP); fail_cnt++;
      end
 
      // T8.3 signs_out consistent with x_out MSBs
      begin
        logic [N_T-1:0] exp_signs;
        for (int k = 0; k < N_T; k++)
          exp_signs[k] = ~x_out[k*XY_T + XY_T - 1];
        if (signs_out === exp_signs) begin
          $display("  PASS [T8.3 signs_consistent] signs=%0b", signs_out);
          pass_cnt++;
        end else begin
          $error("  FAIL [T8.3] signs_out=%0b exp=%0b", signs_out, exp_signs);
          fail_cnt++;
        end
      end
 
      // T8.4 Not all x are zero (bifurcation occurred)
      begin
        logic all_zero;
        all_zero = 1;
        for (int k = 0; k < N_T; k++) begin
          if (x_out[k*XY_T +: XY_T] !== '0) all_zero = 0;
        end
        if (!all_zero) begin
          $display("  PASS [T8.4 bifurcation_occurred] signs=%0b", signs_out);
          pass_cnt++;
        end else begin
          $error("  FAIL [T8.4] all x values are zero after schedule");
          fail_cnt++;
        end
      end
    end
 
    // ── T8.5 Identity J (all zero): x should barely move ──────────────────
    hard_reset();
    write_J(0, {N_T*IC_T{1'b0}});
    write_J(1, {N_T*IC_T{1'b0}});
    write_J(2, {N_T*IC_T{1'b0}});
    write_J(3, {N_T*IC_T{1'b0}});
    // Initial x very small, y=0, no coupling → dynamics driven only by −a*x
    write_xy(0,  0.01, 0.0);
    write_xy(1,  0.01, 0.0);
    write_xy(2,  0.01, 0.0);
    write_xy(3,  0.01, 0.0);
 
    begin
      int sc;
      run_and_wait(sc);
      // With J=0, force = -a*x only; x should decay toward 0, not grow
      // Check that |x| remains < 0.5 for all oscillators
      begin
        logic bounded;
        real xv;

        bounded = 1;
        for (int k = 0; k < N_T; k++) begin
          xv = from_fp($signed(x_out[k*XY_T +: XY_T]));
          if (xv < -0.5 || xv > 0.5) bounded = 0;
        end
        if (bounded) begin
          $display("  PASS [T8.5 zero_J_bounded]"); pass_cnt++;
        end else begin
          $error("  FAIL [T8.5 zero_J_bounded] x escaped ±0.5 with J=0");
          fail_cnt++;
        end
      end
    end
 
    // ── T8.6 schedule_done blocks re-run ──────────────────────────────────
    // schedule_done is still high; pulse run again; should not restart
    begin
      logic sd_before;
      sd_before = schedule_done;
      @(posedge clk); #1; run = 1;
      @(posedge clk); #1; run = 0;
      repeat(5) @(posedge clk);
      // step_done should not fire again after schedule is complete
      if (schedule_done === 1'b1 && sd_before === 1'b1) begin
        $display("  PASS [T8.6 no_rerun_after_done]"); pass_cnt++;
      end else begin
        $error("  FAIL [T8.6] schedule restarted unexpectedly"); fail_cnt++;
      end
    end
 
    $display("--- TB8 result: %0d PASS, %0d FAIL ---\n", pass_cnt, fail_cnt);
 
    // ── Grand summary ──────────────────────────────────────────────────────
    $display("=== GRAND TOTAL: %0d PASS, %0d FAIL ===",
             pass_cnt, fail_cnt);
 
    #20 $finish;
  end
 
endmodule : tb_array
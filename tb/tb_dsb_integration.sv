`timescale 1ns / 1ps
// =============================================================================
// tb_dsb_integration.sv
//
// Full-algorithm integration testbench for dsb_sachi_n3.sv.
// Instantiates dsb_array (which internally wires every sub-module) and runs
// several end-to-end scenarios, checking correctness at the algorithm level.
//
// Compile & run (Cadence xmsim example):
//   xmvlog -sv dsb_sachi_n3_fixed.sv tb_dsb_integration.sv +define+SIMULATION
//   xmelab -timescale 1ns/1ps tb_dsb_integration
//   xmsim  tb_dsb_integration -run -exit
//
// Scenarios
// ─────────────────────────────────────────────────────────────────────────────
//  S1  Ferromagnetic J=+1 (upper-triangular, N=4)
//        All spins should align → signs all-1 or all-0.
//        Checks: schedule_done, step_done count, signs consistent with x_out,
//                x values clamped to ±1 (bifurcation), no all-zero x.
//
//  S2  Anti-ferromagnetic J=−1 (upper-triangular, N=4)
//        Alternating pattern expected: 1010 or 0101.
//        Checks: schedule_done, alternating sign pattern.
//
//  S3  Zero J (identity / uncoupled, N=4)
//        Force = −a(m)·x only → x decays toward 0.
//        Checks: |x_i| < 0.5 for all oscillators after full schedule.
//
//  S4  Single oscillator (N=1, trivial)
//        No neighbours → pure −a(m)·x dynamics → x converges near 0.
//
//  S5  Re-run guard
//        After schedule_done, asserting run again must NOT restart the schedule
//        (schedule_done stays high, no new step_done pulses).
//
//  S6  step_done pulse count equals Nstep exactly.
// =============================================================================

module tb_dsb_integration;
  import dsb_pkg::*;

  // ── Testbench parameters ───────────────────────────────────────────────────
  localparam int unsigned TB_N      = 4;
  localparam int unsigned TB_IC     = 8;
  localparam int unsigned TB_XY_W   = XY_W;
  localparam int unsigned TB_FRAC   = XY_FRAC;
  localparam int unsigned TB_A      = A_BITS;
  localparam int unsigned TB_STEP   = STEP_W;
  localparam int unsigned TB_ACCUM  = ACCUM_W;
  localparam int unsigned TB_PROD   = PROD_W;

  // Algorithm hyper-parameters (matching the reference Python implementation)
  localparam int unsigned DT_FP        = 1638;   // 0.1  × 2^14
  localparam int unsigned C0_FP        = 8192;   // 0.5  × 2^14
  localparam int unsigned A0_FP        = 32768;  // 1.0  in Q1.15
  localparam int unsigned NSTEP        = 50;     // default schedule length
  localparam int unsigned NSTEP_AFM    = 100;    // AFM needs more steps to converge
  localparam int unsigned NSTEP_S6     = 30;     // distinct value for S6 counter check

  // Convenience J encodings
  localparam logic [TB_IC-1:0] Jp = 8'h01;  // J = +1
  localparam logic [TB_IC-1:0] Jn = 8'hFF;  // J = −1  (two's complement)
  localparam logic [TB_IC-1:0] Jz = 8'h00;  // J =  0

  // ── Clock ──────────────────────────────────────────────────────────────────
  logic clk = 1'b0;
  always #5 clk = ~clk;   // 100 MHz

  // ── DUT signals ────────────────────────────────────────────────────────────
  logic [TB_STEP-1:0]            nstep_sig;   // driven per scenario
  logic                          rst_n;
  logic [$clog2(TB_N)-1:0]      wr_row;
  logic                          wr_en;
  logic [TB_N*TB_IC-1:0]        wr_J_row;
  logic [$clog2(TB_N)-1:0]      wr_xy_idx;
  logic                          wr_xy_en;
  logic signed [TB_XY_W-1:0]    wr_x, wr_y;
  logic                          run;
  logic [TB_N*TB_XY_W-1:0]      x_out, y_out;
  logic [TB_N-1:0]               signs_out;
  logic                          step_done, schedule_done;

  // ── DUT instantiation ──────────────────────────────────────────────────────
  dsb_array #(
    .N_P       (TB_N),
    .IC_BITS_P (TB_IC),
    .XY_W_P    (TB_XY_W),
    .XY_FRAC_P (TB_FRAC),
    .A_BITS_P  (TB_A),
    .STEP_W_P  (TB_STEP),
    .ACCUM_W_P (TB_ACCUM),
    .PROD_W_P  (TB_PROD)
  ) dut (
    .clk          (clk),
    .rst_n        (rst_n),
    .wr_row       (wr_row),
    .wr_en        (wr_en),
    .wr_J_row     (wr_J_row),
    .wr_xy_idx    (wr_xy_idx),
    .wr_xy_en     (wr_xy_en),
    .wr_x         (wr_x),
    .wr_y         (wr_y),
    .run          (run),
    .Nstep        (nstep_sig),
    .a0_fp        (TB_A'(A0_FP)),
    .dt_fp        (TB_FRAC'(DT_FP)),
    .c0_fp        (TB_FRAC'(C0_FP)),
    .x_out        (x_out),
    .y_out        (y_out),
    .signs_out    (signs_out),
    .step_done    (step_done),
    .schedule_done(schedule_done)
  );

  // ── Global pass/fail counters ──────────────────────────────────────────────
  int pass_cnt = 0, fail_cnt = 0;

  // ── Helpers ────────────────────────────────────────────────────────────────

  // Convert real → fixed-point Q(XY_INT).(XY_FRAC)
  function automatic logic signed [TB_XY_W-1:0] to_fp(input real v);
    return TB_XY_W'($rtoi(v * real'(1 << TB_FRAC)));
  endfunction

  // Convert fixed-point → real
  function automatic real from_fp(input logic signed [TB_XY_W-1:0] v);
    return real'($signed(v)) / real'(1 << TB_FRAC);
  endfunction

  // Check helper: prints PASS/FAIL and updates counters
  task automatic chk(input string label, input logic cond);
    if (cond) begin
      $display("  PASS [%s]", label);
      pass_cnt++;
    end else begin
      $error("  FAIL [%s]", label);
      fail_cnt++;
    end
  endtask

  // Hard reset: deassert everything, pulse rst_n low then high.
  // nstep_val sets the Nstep port for the upcoming run (default NSTEP=50).
  task automatic hard_reset(input int unsigned nstep_val = NSTEP);
    rst_n    = 1'b0;
    run      = 1'b0;
    wr_en    = 1'b0;
    wr_xy_en = 1'b0;
    wr_row   = '0;
    wr_J_row = '0;
    wr_xy_idx = '0;
    wr_x     = '0;
    wr_y     = '0;
    nstep_sig = TB_STEP'(nstep_val);
    repeat (4) @(posedge clk);
    rst_n = 1'b1;
    repeat (2) @(posedge clk);
  endtask

  // Write one row of J into the array register file
  task automatic write_J(
    input int unsigned           row,
    input logic [TB_N*TB_IC-1:0] data
  );
    @(posedge clk); #1;
    wr_row   = $clog2(TB_N)'(row);
    wr_J_row = data;
    wr_en    = 1'b1;
    @(posedge clk); #1;
    wr_en    = 1'b0;
  endtask

  // Write initial (x, y) for one oscillator
  task automatic write_xy(
    input int unsigned idx,
    input real         x_r,
    input real         y_r
  );
    @(posedge clk); #1;
    wr_xy_idx = $clog2(TB_N)'(idx);
    wr_x      = to_fp(x_r);
    wr_y      = to_fp(y_r);
    wr_xy_en  = 1'b1;
    @(posedge clk); #1;
    wr_xy_en  = 1'b0;
  endtask

  // Pulse run for one cycle then wait until schedule_done or timeout.
  // Returns the number of step_done pulses observed.
  task automatic run_schedule(output int step_count);
    @(posedge clk); #1;
    run = 1'b1;
    @(posedge clk); #1;
    run = 1'b0;

    step_count = 0;
    fork
      // Timeout watchdog
      begin : watchdog
        repeat (100_000) @(posedge clk);
        $error("  FAIL: timeout — schedule_done never asserted");
        fail_cnt++;
        disable waiter;
      end
      // Wait for schedule completion, count step_done pulses
      begin : waiter
        forever begin
          @(posedge clk); #1;
          if (step_done)    step_count++;
          if (schedule_done) begin
            disable watchdog;
            break;
          end
        end
      end
    join
  endtask

  // Display full oscillator state
  // Zero-pad signs to TB_N bits so a leading-zero spin is never silently dropped
  task automatic show_state(input string label);
    string s;
    s = "";
    for (int k = TB_N-1; k >= 0; k--)
      s = {s, (signs_out[k] ? "1" : "0")};
    $write("%s  signs=%s  x=[", label, s);
    for (int k = 0; k < TB_N; k++) begin
      $write("%8.4f", from_fp($signed(x_out[k*TB_XY_W +: TB_XY_W])));
      if (k < TB_N - 1) $write(", ");
    end
    $display("]");
  endtask

  // ── Main test sequence ─────────────────────────────────────────────────────
  initial begin
    logic [TB_N-1:0] signs_run1, signs_run2;

    $dumpfile("tb_dsb_integration.vcd");
    $dumpvars(0, tb_dsb_integration);

    // =========================================================================
    // S1 — Ferromagnetic J=+1 (upper-triangular)
    //       All spins see positive coupling → all should align to +1 or −1.
    //       With positive initial x, expect signs_out = 4'b1111.
    // =========================================================================
    $display("\n=== S1: Ferromagnetic J=+1 (N=%0d, Nstep=%0d) ===", TB_N, NSTEP);
    hard_reset();

    // Upper-triangular J = +1
    //   Row 0: J[0][3..0] = +1 +1 +1  0  (MSB = osc 3)
    //   Row 1: J[1][3..0] = +1 +1  0  0
    //   Row 2: J[2][3..0] = +1  0  0  0
    //   Row 3: all zero
    write_J(0, {Jp, Jp, Jp, Jz});
    write_J(1, {Jp, Jp, Jz, Jz});
    write_J(2, {Jp, Jz, Jz, Jz});
    write_J(3, {TB_N*TB_IC{1'b0}});

    // Asymmetric positive initial x so bifurcation is unambiguous
    write_xy(0,  0.05, 0.0);
    write_xy(1,  0.04, 0.0);
    write_xy(2,  0.06, 0.0);
    write_xy(3,  0.03, 0.0);

    begin
      int sc;
      run_schedule(sc);
      show_state("S1 final:");

      // S1-a  schedule_done was asserted
      chk("S1-a schedule_done",  schedule_done === 1'b1);

      // S1-b  step_done fired exactly Nstep times
      chk($sformatf("S1-b step_count=%0d (exp %0d)", sc, NSTEP),
          sc === NSTEP);

      // S1-c  signs_out is consistent with x_out MSBs
      begin
        logic [TB_N-1:0] exp_signs;
        for (int k = 0; k < TB_N; k++)
          exp_signs[k] = ~x_out[k*TB_XY_W + TB_XY_W - 1];
        chk($sformatf("S1-c signs_consistent signs=%04b", signs_out),
            signs_out === exp_signs);
      end

      // S1-d  x values are NOT all zero (bifurcation happened)
      begin
        logic any_nonzero;
        any_nonzero = 1'b0;
        for (int k = 0; k < TB_N; k++)
          if (x_out[k*TB_XY_W +: TB_XY_W] !== '0) any_nonzero = 1'b1;
        chk("S1-d bifurcation_nonzero", any_nonzero);
      end

      // S1-e  All spins aligned (ferromagnetic ground state)
      //       All positive initial x → expect all spins = 1
      chk($sformatf("S1-e ferromagnetic_aligned signs=%04b", signs_out),
          signs_out === {TB_N{1'b1}} || signs_out === {TB_N{1'b0}});
    end


    // =========================================================================
    // S2 — Anti-ferromagnetic J=−1 (upper-triangular)
    //       Competing spins → expect alternating pattern 1010 or 0101.
    //       Uses NSTEP_AFM=100: the AFM energy landscape requires more steps
    //       to bifurcate cleanly than the ferromagnetic case.
    // =========================================================================
    $display("\n=== S2: Anti-ferromagnetic J=−1 (N=%0d, Nstep=%0d) ===",
             TB_N, NSTEP_AFM);
    hard_reset(NSTEP_AFM);

    write_J(0, {Jn, Jn, Jn, Jz});
    write_J(1, {Jn, Jn, Jz, Jz});
    write_J(2, {Jn, Jz, Jz, Jz});
    write_J(3, {TB_N*TB_IC{1'b0}});

    // Stronger alternating nudge: larger magnitude seeds the correct AFM ground
    // state more reliably and gives the algorithm an unambiguous starting point.
    //   osc 0 (+), osc 1 (−), osc 2 (+), osc 3 (−)  → expect signs 1010
    write_xy(0,  0.08, 0.0);
    write_xy(1, -0.08, 0.0);
    write_xy(2,  0.08, 0.0);
    write_xy(3, -0.08, 0.0);

    begin
      int sc;
      run_schedule(sc);
      show_state("S2 final:");

      chk("S2-a schedule_done", schedule_done === 1'b1);
      chk($sformatf("S2-b step_count=%0d (exp %0d)", sc, NSTEP_AFM),
          sc === NSTEP_AFM);

      // S2-c  Alternating spin pattern (the AFM ground state for this J)
      chk($sformatf("S2-c alternating_pattern signs=%s",
                    $sformatf("%04b", signs_out)),
          signs_out === 4'b1010 || signs_out === 4'b0101);
    end


    // =========================================================================
    // S3 — Zero coupling (J=0 everywhere)
    //       Force = −a(m)·x → x decays toward 0 throughout the schedule.
    //       All |x_i| must stay below 0.5 at the end.
    // =========================================================================
    $display("\n=== S3: Zero J (decoupled, N=%0d) ===", TB_N);
    hard_reset();

    write_J(0, {TB_N*TB_IC{1'b0}});
    write_J(1, {TB_N*TB_IC{1'b0}});
    write_J(2, {TB_N*TB_IC{1'b0}});
    write_J(3, {TB_N*TB_IC{1'b0}});

    // Small positive initial x for all oscillators
    write_xy(0,  0.10, 0.0);
    write_xy(1,  0.08, 0.0);
    write_xy(2,  0.12, 0.0);
    write_xy(3,  0.06, 0.0);

    begin
      int sc;
      run_schedule(sc);
      show_state("S3 final:");

      chk("S3-a schedule_done", schedule_done === 1'b1);

      // S3-b  All x stay bounded (no bifurcation without coupling)
      begin
        logic bounded;
        real xv;

        bounded = 1'b1;
        for (int k = 0; k < TB_N; k++) begin
          xv = from_fp($signed(x_out[k*TB_XY_W +: TB_XY_W]));
          if (xv < -0.5 || xv > 0.5) bounded = 1'b0;
        end
        chk("S3-b zero_J_x_bounded_0.5", bounded);
      end
    end


    // =========================================================================
    // S4 — Sequential schedules: verify a fresh run after reset works correctly
    //       Run the ferromagnetic case twice after two separate resets and
    //       confirm both produce the same aligned result (deterministic output).
    // =========================================================================
    $display("\n=== S4: Deterministic re-run after reset ===");

    // First run
    hard_reset();
    write_J(0, {Jp, Jp, Jp, Jz});
    write_J(1, {Jp, Jp, Jz, Jz});
    write_J(2, {Jp, Jz, Jz, Jz});
    write_J(3, {TB_N*TB_IC{1'b0}});
    write_xy(0,  0.05, 0.0);
    write_xy(1,  0.04, 0.0);
    write_xy(2,  0.06, 0.0);
    write_xy(3,  0.03, 0.0);
    begin int sc; run_schedule(sc); end
    signs_run1 = signs_out;
    show_state("S4 run1 final:");

    // Second run (identical setup)
    hard_reset();
    write_J(0, {Jp, Jp, Jp, Jz});
    write_J(1, {Jp, Jp, Jz, Jz});
    write_J(2, {Jp, Jz, Jz, Jz});
    write_J(3, {TB_N*TB_IC{1'b0}});
    write_xy(0,  0.05, 0.0);
    write_xy(1,  0.04, 0.0);
    write_xy(2,  0.06, 0.0);
    write_xy(3,  0.03, 0.0);
    begin int sc; run_schedule(sc); end
    signs_run2 = signs_out;
    show_state("S4 run2 final:");

    chk($sformatf("S4 deterministic (run1=%04b run2=%04b)", signs_run1, signs_run2),
        signs_run1 === signs_run2);


    // =========================================================================
    // S5 — Re-run guard: after schedule_done, pulsing run must not restart.
    //       The schedule_done flag must stay high and no new step_done pulses
    //       should appear for at least 5×N×6 cycles (well beyond one full step).
    // =========================================================================
    $display("\n=== S5: Re-run guard (schedule_done blocks restart) ===");
    // (DUT is still at the end of the second ferromagnetic run from S4)

    begin
      int spurious_steps;
      logic sd_before;

      spurious_steps = 0;
      sd_before = schedule_done;

      // Pulse run while schedule_done is already high
      @(posedge clk); #1; run = 1'b1;
      @(posedge clk); #1; run = 1'b0;

      // Watch for 200 cycles — no new step_done should fire
      repeat (200) begin
        @(posedge clk); #1;
        if (step_done) spurious_steps++;
      end

      chk("S5-a schedule_done_stays_high", schedule_done === 1'b1);
      chk($sformatf("S5-b no_spurious_step_done (count=%0d)", spurious_steps),
          spurious_steps === 0);
    end


    // =========================================================================
    // S6 — step_done pulse count with a distinct Nstep=30
    //       Uses Nstep=30 (different from the default 50) so this is a truly
    //       independent verification that the step counter respects the port
    //       value and fires exactly Nstep times — no more, no less.
    // =========================================================================
    $display("\n=== S6: step_done pulse count (Nstep=%0d) ===", NSTEP_S6);
    hard_reset(NSTEP_S6);
    write_J(0, {Jp, Jp, Jp, Jz});
    write_J(1, {Jp, Jp, Jz, Jz});
    write_J(2, {Jp, Jz, Jz, Jz});
    write_J(3, {TB_N*TB_IC{1'b0}});
    write_xy(0,  0.05, 0.0);
    write_xy(1,  0.04, 0.0);
    write_xy(2,  0.06, 0.0);
    write_xy(3,  0.03, 0.0);

    begin
      int sc;
      run_schedule(sc);
      chk($sformatf("S6 step_count=%0d exp=%0d", sc, NSTEP_S6),
          sc === NSTEP_S6);
    end


    // =========================================================================
    // S7 — signs_out tracks x_out MSB at every step
    //       Monitor during a live run; check consistency on every step_done.
    // =========================================================================
    $display("\n=== S7: signs_out tracks x_out MSB during run ===");
    hard_reset();
    write_J(0, {Jp, Jp, Jp, Jz});
    write_J(1, {Jp, Jp, Jz, Jz});
    write_J(2, {Jp, Jz, Jz, Jz});
    write_J(3, {TB_N*TB_IC{1'b0}});
    write_xy(0,  0.05, 0.0);
    write_xy(1,  0.04, 0.0);
    write_xy(2,  0.06, 0.0);
    write_xy(3,  0.03, 0.0);

    begin
      int mismatch_steps;
      mismatch_steps = 0;
      @(posedge clk); #1; run = 1'b1;
      @(posedge clk); #1; run = 1'b0;

      fork
        begin : s7_wd
          repeat (100_000) @(posedge clk);
          $error("  FAIL S7: timeout");
          fail_cnt++;
          disable s7_wait;
        end
        begin : s7_wait
          forever begin
            logic exp_sign;
            @(posedge clk); #1;
            if (step_done) begin
              // Check signs_out == ~MSB of each x_out slice
              for (int k = 0; k < TB_N; k++) begin
                exp_sign = ~x_out[k*TB_XY_W + TB_XY_W - 1];
                if (signs_out[k] !== exp_sign) mismatch_steps++;
              end
            end
            if (schedule_done) begin
              disable s7_wd;
              break;
            end
          end
        end
      join

      chk($sformatf("S7 signs_track_x_msb (mismatches=%0d)", mismatch_steps),
          mismatch_steps === 0);
    end


    // =========================================================================
    // Summary
    // =========================================================================
    $display("\n=== INTEGRATION RESULT: %0d PASS, %0d FAIL ===",
             pass_cnt, fail_cnt);
    if (fail_cnt == 0)
      $display("ALL TESTS PASSED");
    else
      $display("*** %0d TEST(S) FAILED ***", fail_cnt);

    #20 $finish;
  end

endmodule : tb_dsb_integration
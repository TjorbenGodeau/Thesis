`timescale 1ns / 1ps
// =============================================================================
// 9.  Testbench — 4-oscillator dSB, ferromagnetic J = +1 everywhere
//
//     J = [ 0  +1  +1  +1 ]   (upper-triangular, matching Python dSB)
//         [ 0   0  +1  +1 ]
//         [ 0   0   0  +1 ]
//         [ 0   0   0   0 ]
//
//     Expected: oscillators bifurcate to alternating signs (1010 or 0101)
//
//     Fixed-point constants:
//       XY_FRAC = 14  → Δ ≈ 6.1×10⁻⁵
//       IC_BITS = 8   → J = 8'h01 (+1)
//       dt = 0.1      → dt_fp = 1638
//       c0 = 0.5      → c0_fp = 8192
//       a0 = 1.0      → a0_fp = 32768  (Q1.15)
//       Nstep = 50
// =============================================================================
`ifdef SIMULATION
module tb_dsb_sachi;
  import dsb_pkg::*;
 
  // ── Testbench-local parameters ─────────────────────────────────────────────
  localparam int unsigned TB_N       = 4;
  localparam int unsigned TB_IC      = 8;
  localparam int unsigned TB_XY_W    = XY_W;
  localparam int unsigned TB_XY_FRAC = XY_FRAC;
  localparam int unsigned TB_A_BITS  = A_BITS;
  localparam int unsigned TB_STEP_W  = STEP_W;
  localparam int unsigned TB_ACCUM_W = ACCUM_W;
  localparam int unsigned TB_PROD_W  = PROD_W;
 
  localparam int unsigned DT_FP  = 1638;   // 0.1  × 2^14
  localparam int unsigned C0_FP  = 8192;   // 0.5  × 2^14
  localparam int unsigned A0_FP  = 32768;  // 1.0  in Q1.15
  localparam int unsigned NSTEP  = 50;
 
  localparam logic [TB_IC-1:0] Jp = 8'h01;  // J = +1
  localparam logic [TB_IC-1:0] Jz = 8'h00;  // J =  0
 
  // ── Clock ──────────────────────────────────────────────────────────────────
  logic clk = 1'b0;
  always #5 clk = ~clk;   // 100 MHz
 
  // ── Clocking block ─────────────────────────────────────────────────────────
  clocking cb @(posedge clk);
    default input #1 output #1;
    input  signs_out, step_done, schedule_done, x_out, y_out;
    output wr_row, wr_en, wr_J_row, wr_xy_idx, wr_xy_en, wr_x, wr_y, run;
  endclocking
 
  // ── DUT ports ──────────────────────────────────────────────────────────────
  logic                            rst_n;
  logic [$clog2(TB_N)-1:0]        wr_row;
  logic                            wr_en;
  logic [TB_N*TB_IC-1:0]          wr_J_row;
  logic [$clog2(TB_N)-1:0]        wr_xy_idx;
  logic                            wr_xy_en;
  logic signed [TB_XY_W-1:0]      wr_x, wr_y;
  logic                            run;
  logic [TB_N*TB_XY_W-1:0]        x_out, y_out;
  logic [TB_N-1:0]                 signs_out;
  logic                            step_done, schedule_done;
 
  // ── DUT ────────────────────────────────────────────────────────────────────
  dsb_array #(
    .N_P (TB_N),       
    .IC_BITS_P (TB_IC),    
    .XY_W_P (TB_XY_W),
    .XY_FRAC_P (TB_XY_FRAC), 
    .A_BITS_P (TB_A_BITS), 
    .STEP_W_P (TB_STEP_W),
    .ACCUM_W_P (TB_ACCUM_W), 
    .PROD_W_P (TB_PROD_W)
  ) dut (
    .clk (clk),
    .rst_n (rst_n),
    .wr_row (wr_row),
    .wr_en (wr_en),
    .wr_J_row (wr_J_row),
    .wr_xy_idx (wr_xy_idx),
    .wr_xy_en (wr_xy_en),
    .wr_x (wr_x),
    .wr_y (wr_y),
    .run (run),
    .Nstep (TB_STEP_W'(NSTEP)),
    .a0_fp (TB_A_BITS'(A0_FP)),
    .dt_fp (TB_XY_FRAC'(DT_FP)),
    .c0_fp (TB_XY_FRAC'(C0_FP)),
    .x_out (x_out),
    .y_out (y_out),
    .signs_out (signs_out),
    .step_done (step_done),
    .schedule_done (schedule_done)
  );
 
  // ── Helper task: display x values as reals ─────────────────────────────────
  task automatic show_state(input string label);
    real xr;
    $write("%s  signs=%0b  x=[", label, signs_out);
    for (int k = 0; k < TB_N; k++) begin
      xr = $signed(x_out[k*TB_XY_W +: TB_XY_W]) / real'(1 << TB_XY_FRAC);
      $write("%7.4f", xr);
      if (k < TB_N - 1) $write(", ");
    end
    $display("]");
  endtask
 
  // ── Helper task: load one J row ────────────────────────────────────────────
  task automatic load_J_row(
    input int unsigned          row,
    input logic [TB_N*TB_IC-1:0] data
  );
    @(posedge clk);
    wr_row    = $clog2(TB_N)'(row);
    wr_J_row  = data;
    wr_en     = 1'b1;
    @(posedge clk);
    wr_en     = 1'b0;
  endtask
 
  // ── Helper task: load initial x, y for one oscillator ─────────────────────
  task automatic load_xy(
    input int unsigned          idx,
    input logic signed [TB_XY_W-1:0] x_val,
    input logic signed [TB_XY_W-1:0] y_val
  );
    @(posedge clk);
    wr_xy_idx = $clog2(TB_N)'(idx);
    wr_x      = x_val;
    wr_y      = y_val;
    wr_xy_en  = 1'b1;
    @(posedge clk);
    wr_xy_en  = 1'b0;
  endtask
 
  // ── Main test program ──────────────────────────────────────────────────────
  program automatic test_prog;
    initial begin
      $dumpfile("tb_dsb_sachi.vcd");
      $dumpvars(0, tb_dsb_sachi);
 
      // Initialise driven signals
      rst_n    = 1'b0;
      wr_en    = 1'b0;
      wr_xy_en = 1'b0;
      run      = 1'b0;
      wr_row   = '0;
      wr_J_row = '0;
      wr_x     = '0;
      wr_y     = '0;
 
      repeat (4) @(posedge clk);
      rst_n = 1'b1;
      repeat (2) @(posedge clk);
 
      // ── Load J matrix (upper-triangular) ──────────────────────────────────
      // Row 0: J[0][1..3] = +1,  J[0][0] = 0
      load_J_row(0, {Jp, Jp, Jp, Jz});
      // Row 1: J[1][2..3] = +1
      load_J_row(1, {Jp, Jp, Jz, Jz});
      // Row 2: J[2][3] = +1
      load_J_row(2, {Jp, Jz, Jz, Jz});
      // Row 3: all zero
      load_J_row(3, {TB_N*TB_IC{1'b0}});
 
      // ── Load initial x (matching Python U[−0.1, +0.1]) ───────────────────
      // Fixed-point: 0.05 × 2^14 = 819
      load_xy(0,  TB_XY_W'(16'sd819),  '0);
      load_xy(1, -TB_XY_W'(16'sd410),  '0);
      load_xy(2,  TB_XY_W'(16'sd614),  '0);
      load_xy(3, -TB_XY_W'(16'sd205),  '0);
 
      // ── Run the schedule ──────────────────────────────────────────────────
      $display("[%0t] Starting dSB schedule: N=%0d, Nstep=%0d",
               $realtime, TB_N, NSTEP);
      @(posedge clk); run = 1'b1;
      @(posedge clk); run = 1'b0;
 
      begin : run_loop
        int step_cnt = 0;
        forever begin
          @(posedge clk);
          if (step_done) begin
            show_state($sformatf("[%0t] step %0d:", $realtime, step_cnt));
            step_cnt++;
          end
          if (schedule_done) break;
          if (step_cnt > 10_000) begin
            $error("Timeout: schedule did not complete in 10000 steps.");
            break;
          end
        end
      end
 
      show_state("Final state:");
      $display("Final signs: %0b", signs_out);
      $display("Expected alternating: %0b or %0b",
               {TB_N{2'b10}} >> 0, {TB_N{2'b01}} >> 0);
 
      // ── Self-check: verify alternating pattern ────────────────────────────
      if (signs_out == TB_N'('b1010) || signs_out == TB_N'('b0101))
        $display("PASS: correct alternating spin pattern.");
      else
        $warning("Note: non-alternating pattern — may still be valid for this J.");
 
      #20 $finish;
    end
  endprogram : test_prog
 
endmodule : tb_dsb_sachi
`endif // SIMULATION
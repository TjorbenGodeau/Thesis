// tb_update_unit.sv
`timescale 1ns/1ps

module tb_update_unit;

    import dsb_pkg::*;

    // ── Test parameters ────────────────────────────────────────────────────
    localparam int XY_W_P    = XY_W;
    localparam int XY_FRAC_P = XY_FRAC;
    localparam int A_BITS_P  = A_BITS;
    localparam int ACCUM_W_P = ACCUM_W;
    localparam int PROD_W_P  = PROD_W;

    // Fixed‑point constants
    localparam logic signed [XY_W_P-1:0] ONE_FP_LOCAL = XY_W_P'(1 << XY_FRAC_P);

    // Clock
    logic clk = 0;
    always #5 clk = ~clk;

    // DUT signals
    logic                          start;
    logic signed [XY_W_P-1:0]      x_i;
    logic signed [XY_W_P-1:0]      y_i;
    logic signed [ACCUM_W_P-1:0]   Jx_i;
    logic [A_BITS_P-1:0]           a_m;
    logic [XY_FRAC_P-1:0]          dt_fp;
    logic [XY_FRAC_P-1:0]          c0_fp;
    logic signed [XY_W_P-1:0]      x_i_new;
    logic signed [XY_W_P-1:0]      y_i_new;
    logic                          done;

    // DUT
    update_unit u_dut (
        .clk     (clk),
        .start   (start),
        .x_i     (x_i),
        .y_i     (y_i),
        .Jx_i    (Jx_i),
        .a_m     (a_m),
        .dt_fp   (dt_fp),
        .c0_fp   (c0_fp),
        .x_i_new (x_i_new),
        .y_i_new (y_i_new),
        .done    (done)
    );

    // ── Reference model (same arithmetic as RTL) ──────────────────────────
    function automatic void ref_update(
        input  logic signed [XY_W_P-1:0]      x_i,
        input  logic signed [XY_W_P-1:0]      y_i,
        input  logic signed [ACCUM_W_P-1:0]   Jx_i,
        input  logic [A_BITS_P-1:0]           a_m,
        input  logic [XY_FRAC_P-1:0]          dt_fp,
        input  logic [XY_FRAC_P-1:0]          c0_fp,
        output logic signed [XY_W_P-1:0]      x_ref,
        output logic signed [XY_W_P-1:0]      y_ref
    );
        // Stage A: a_m * x_i
        logic signed [PROD_W_P-1:0]   ax_full;
        logic signed [XY_W_P-1:0]     ax_shifted;
        ax_full = $signed({{(PROD_W_P-XY_W_P){x_i[XY_W_P-1]}}, x_i})
                  * $signed({1'b0, a_m});
        ax_shifted = XY_W_P'(ax_full >>> (A_BITS_P - 1));

        // Stage B: c0 * Jx_i
        logic signed [ACCUM_W_P+XY_FRAC_P-1:0] c0jx_full;
        c0jx_full = Jx_i * $signed({1'b0, c0_fp});

        // Stage C: Δy = dt * (-a*x + c0*Jx)
        logic signed [XY_W_P-1:0]   neg_ax;
        logic signed [ACCUM_W_P+XY_FRAC_P-1:0] force_wide;
        logic signed [ACCUM_W_P+2*XY_FRAC_P-1:0] dy_full;
        logic signed [XY_W_P-1:0]   delta_y;
        neg_ax = -ax_shifted;
        force_wide = $signed({{(ACCUM_W_P+XY_FRAC_P-XY_W_P){neg_ax[XY_W_P-1]}}, neg_ax})
                     + $signed(c0jx_full);
        dy_full = $signed(force_wide) * $signed({1'b0, dt_fp});
        delta_y = XY_W_P'(dy_full >>> XY_FRAC_P);

        // Stage D: y_pre = y_i + delta_y
        logic signed [XY_W_P-1:0] y_pre;
        y_pre = y_i + delta_y;

        // Stage E: x_pre = x_i + dt * y_pre
        logic signed [2*XY_W_P-1:0] dx_full;
        logic signed [XY_W_P-1:0]   x_pre;
        dx_full = $signed(y_pre) * $signed({1'b0, dt_fp});
        x_pre = x_i + XY_W_P'(dx_full >>> XY_FRAC_P);

        // Stage F: wall clamp
        logic wall_hit_pos, wall_hit_neg, wall_hit;
        wall_hit_pos = (x_pre >  $signed(ONE_FP_LOCAL));
        wall_hit_neg = (x_pre < -$signed(ONE_FP_LOCAL));
        wall_hit = wall_hit_pos | wall_hit_neg;

        if (wall_hit) begin
            x_ref = wall_hit_pos ?  ONE_FP_LOCAL : -ONE_FP_LOCAL;
            y_ref = '0;
        end else begin
            x_ref = x_pre;
            y_ref = y_pre;
        end
    endfunction

    // ── Test control ──────────────────────────────────────────────────────
    int error_count = 0;

    task check_equal(input string msg, input logic signed [XY_W_P-1:0] a,
                     input logic signed [XY_W_P-1:0] b);
        if (a !== b) begin
            $error("%s: expected %0d, got %0d", msg, b, a);
            error_count++;
        end else begin
            $display("  %s: OK", msg);
        end
    endtask

    task test_case(
        input string label,
        input logic signed [XY_W_P-1:0] x_i_in,
        input logic signed [XY_W_P-1:0] y_i_in,
        input logic signed [ACCUM_W_P-1:0] Jx_i_in,
        input logic [A_BITS_P-1:0] a_m_in,
        input logic [XY_FRAC_P-1:0] dt_fp_in,
        input logic [XY_FRAC_P-1:0] c0_fp_in
    );
        logic signed [XY_W_P-1:0] x_exp, y_exp;
        // Compute reference
        ref_update(x_i_in, y_i_in, Jx_i_in, a_m_in, dt_fp_in, c0_fp_in, x_exp, y_exp);

        // Apply inputs
        x_i  = x_i_in;
        y_i  = y_i_in;
        Jx_i = Jx_i_in;
        a_m  = a_m_in;
        dt_fp = dt_fp_in;
        c0_fp = c0_fp_in;

        // Start the update
        @(posedge clk);
        start = 1;
        @(posedge clk);
        start = 0;

        // done should be high next cycle
        @(posedge clk);
        if (!done) $error("%s: done not asserted", label);
        @(posedge clk);  // extra cycle for outputs to settle

        // Check outputs
        check_equal($sformatf("%s x_i_new", label), x_i_new, x_exp);
        check_equal($sformatf("%s y_i_new", label), y_i_new, y_exp);
    endtask

    // ── Test sequence ──────────────────────────────────────────────────────
    initial begin
        start = 0;
        x_i = 0; y_i = 0; Jx_i = 0; a_m = 0; dt_fp = 0; c0_fp = 0;

        repeat (2) @(posedge clk);
        $display("=== Starting update_unit test ===");

        // ── Test 1: Zero inputs → no change ──────────────────────────────────
        $display("Test 1: All zero");
        test_case("Test1", 0, 0, 0, 0, 0, 0);

        // ── Test 2: Pure momentum (a=1, dt=0.01, no Jx) ─────────────────────
        // a=1 (Q0.14: 16384), dt=0.01 (164), x=100, y=0, Jx=0.
        // Expected: dx = 0 (since y=0), dy = -a*x*dt = -100*0.01 = -1, y_pre = -1, x_pre = 100.
        $display("Test 2: Pure momentum (a=1, x=100)");
        test_case("Test2", 100, 0, 0, 16384, 164, 0);

        // ── Test 3: Jx coupling only (Jx=50, c0=0.5, dt=0.01, x=0, y=0) ─────
        // force = c0*Jx = 0.5*50 = 25; dy = 25*0.01 = 0.25 → in Q0.14: 0.25*16384=4096
        // y_pre = 4096 (in Q0.14), dx = y_pre*dt = 4096*0.01 = 40.96 ~ 40 (truncated)
        // x_pre = 40.
        $display("Test 3: Jx coupling (Jx=50, c0=0.5)");
        test_case("Test3", 0, 0, 50, 0, 164, 8192);  // 0.5*16384=8192

        // ── Test 4: Combined (a=0.5, Jx=-20, c0=0.3, dt=0.01, x=10, y=5) ────
        // Use reference function to compute expected manually? We rely on ref.
        $display("Test 4: Combined case");
        test_case("Test4", 10, 5, -20, 8192, 164, 4915); // 0.3*16384≈4915

        // ── Test 5: Wall hit (x_pre > 1) ─────────────────────────────────────
        // Large Jx and c0 to push x over 1.
        $display("Test 5: Wall hit positive");
        test_case("Test5", 0, 0, 200, 0, 164, 16384); // c0=1, dt=0.01, Jx=200 → force=200, dy=2, x_pre≈200*dt? Actually x_pre = dt*y_pre, y_pre=2, so x_pre≈0.02, not enough. Need larger dt.
        // Let's use dt=0.1 (1638 in Q0.14), Jx=500, c0=1 → force=500, dy=50, x_pre=5 → hit +1.
        $display("Test 5b: Wall hit positive with dt=0.1");
        test_case("Test5b", 0, 0, 500, 0, 1638, 16384); // dt_fp = 0.1*16384=1638

        // ── Test 6: Wall hit negative ─────────────────────────────────────
        test_case("Test6", 0, 0, -500, 0, 1638, 16384);

        // ── Test 7: a_m zero → no momentum, but Jx coupling ──────────────
        test_case("Test7", 50, 10, 30, 0, 164, 8192);

        // ── Test 8: Large a_m with negative x ─────────────────────────────
        test_case("Test8", -50, 5, 10, 16384, 164, 0);

        // ── Test 9: dt=0 → no change ──────────────────────────────────────
        test_case("Test9", 20, -10, 100, 16384, 0, 8192);

        // ── Summary ──────────────────────────────────────────────────────────
        if (error_count == 0) begin
            $display("=== All tests PASSED ===");
        end else begin
            $error("=== %0d tests FAILED ===", error_count);
            $finish(1);
        end
        $finish;
    end

    // VCD dump
    initial begin
        $dumpfile("tb_update_unit.vcd");
        $dumpvars(0, tb_update_unit);
    end

endmodule
// tb_dotprod_phase2_sparse.sv – Self-checking testbench for the corrected dotprod
`timescale 1ns/1ps

module tb_dotprod_phase2_sparse;

    import dsb_pkg::*;
    import dsb_sparse_pkg::*;

    localparam int K_MAX_P   = 4;
    localparam int IC_BITS_P = 4;
    localparam int ACCUM_W_P = 16;
    localparam int K_CNT_W   = $clog2(K_MAX_P + 1);

    logic clk = 0;
    always #5 clk = ~clk;

    logic                         start;
    logic                         done;
    logic                         accumulate;
    logic                         last_chunk;
    logic [K_CNT_W-1:0]           valid_count;
    logic                         sign_xi;
    logic [K_MAX_P*IC_BITS_P-1:0] xnor_J;
    logic [K_MAX_P-1:0]           sign_eq;
    logic signed [ACCUM_W_P-1:0]  Jx_i;

    dotprod_phase2_sparse #(
        .K_MAX_P   (K_MAX_P),
        .IC_BITS_P (IC_BITS_P),
        .ACCUM_W_P (ACCUM_W_P)
    ) u_dut (
        .clk         (clk),
        .start       (start),
        .done        (done),
        .accumulate  (accumulate),
        .last_chunk  (last_chunk),
        .valid_count (valid_count),
        .sign_xi     (sign_xi),
        .xnor_J      (xnor_J),
        .sign_eq     (sign_eq),
        .Jx_i        (Jx_i)
    );

    int error_count = 0;

    // Helper to pack a single slot's xnor_J into the vector
    function automatic logic [K_MAX_P*IC_BITS_P-1:0] pack_xnor_J(
        input int slot, input logic [IC_BITS_P-1:0] value
    );
        return value << (slot * IC_BITS_P);
    endfunction

    // Reference model – implements the same per‑term logic as the corrected RTL
    function automatic logic signed [ACCUM_W_P-1:0] ref_dotprod(
        input logic [ACCUM_W_P-1:0] prev_accum,
        input logic                  accumulate,
        input logic [K_CNT_W-1:0]    valid_count,
        input logic                  sign_xi,
        input logic [K_MAX_P*IC_BITS_P-1:0] xnor_J,
        input logic [K_MAX_P-1:0]    sign_eq
    );
        logic signed [ACCUM_W_P-1:0] chunk_sum;
        chunk_sum = '0;
        for (int k = 0; k < K_MAX_P; k++) begin
            if (k < int'(valid_count)) begin
                logic [IC_BITS_P-1:0] xnor_word;
                logic signed [IC_BITS_P-1:0] J_val;
                xnor_word = xnor_J[k*IC_BITS_P +: IC_BITS_P];
                J_val = sign_xi ? signed'(xnor_word) : signed'(~xnor_word);
                chunk_sum += ACCUM_W_P'(sign_eq[k] ? J_val : -J_val);
            end
        end
        return (accumulate ? prev_accum : '0) + chunk_sum;
    endfunction

    // Run a single test case
    task automatic run_test(
        input string label,
        input logic  accumulate_in,
        input logic  last_chunk_in,
        input logic [K_CNT_W-1:0] valid_count_in,
        input logic  sign_xi_in,
        input logic [K_MAX_P*IC_BITS_P-1:0] xnor_J_in,
        input logic [K_MAX_P-1:0] sign_eq_in,
        input logic signed [ACCUM_W_P-1:0] expected
    );
        accumulate  = accumulate_in;
        last_chunk  = last_chunk_in;
        valid_count = valid_count_in;
        sign_xi     = sign_xi_in;
        xnor_J      = xnor_J_in;
        sign_eq     = sign_eq_in;

        @(posedge clk);
        start = 1;
        @(posedge clk);
        start = 0;

        // The DUT asserts done on this same edge; check immediately.
        if (!done) $error("%s: done not asserted", label);

        if (Jx_i !== expected) begin
            $error("%s: expected %0d, got %0d", label, expected, Jx_i);
            error_count++;
        end else begin
            $display("%s: OK (Jx_i=%0d)", label, Jx_i);
        end
    endtask

    // Test sequence
    initial begin
        // ---- all local declarations at top ----
        logic [K_MAX_P*IC_BITS_P-1:0] xnor1, xnor_chunk0, xnor_chunk1;
        logic [K_MAX_P*IC_BITS_P-1:0] xnor_c0_neg, xnor_c1_neg, xnor_partial;

        start       = 0;
        accumulate  = 0;
        last_chunk  = 0;
        valid_count = 0;
        sign_xi     = 0;
        xnor_J      = 0;
        sign_eq     = 0;

        repeat (2) @(posedge clk);
        $display("=== Starting dotprod_phase2_sparse test (reference model) ===");

        // --------------------------------------------------------------------
        // Test 1: sign_xi=+1, all sign_eq=1, J=[5,-3,2,1] -> sum=5
        // --------------------------------------------------------------------
        $display("Test 1: sign_xi=+1, all sign_eq=1");
        xnor1 = 0;
        xnor1 |= pack_xnor_J(0, 4'sd5);
        xnor1 |= pack_xnor_J(1, -4'sd3);
        xnor1 |= pack_xnor_J(2, 4'sd2);
        xnor1 |= pack_xnor_J(3, 4'sd1);
        run_test("Test1", 0, 1, 4, 1, xnor1, 4'b1111, 5);

        // --------------------------------------------------------------------
        // Test 2: sign_xi=-1, all sign_eq=0, same J.
        // With sign_xi=-1, sign_eq=0 means sign_j = +1 (opposite).
        // So term = J * (+1) * (-1) = -J. Sum = -5 +3 -2 -1 = -5.
        // --------------------------------------------------------------------
        $display("Test 2: sign_xi=-1, all sign_eq=0");
        run_test("Test2", 0, 1, 4, 0, xnor1, 4'b0000, -5);

        // --------------------------------------------------------------------
        // Test 3: Two chunks, sign_xi=+1, J=[1,2,3,4] then [5,6] -> total 21
        // --------------------------------------------------------------------
        $display("Test 3: Two chunks, sign_xi=+1");
        xnor_chunk0 = 0;
        xnor_chunk1 = 0;
        xnor_chunk0 |= pack_xnor_J(0, 4'sd1) | pack_xnor_J(1, 4'sd2) |
                       pack_xnor_J(2, 4'sd3) | pack_xnor_J(3, 4'sd4);
        xnor_chunk1 |= pack_xnor_J(0, 4'sd5) | pack_xnor_J(1, 4'sd6);
        run_test("Test3a chunk0", 0, 0, 4, 1, xnor_chunk0, 4'b1111, 10);
        run_test("Test3b chunk1", 1, 1, 2, 1, xnor_chunk1, 2'b11, 21);

        // --------------------------------------------------------------------
        // Test 4: Two chunks, sign_xi=-1, all sign_eq=0.
        // For sign_xi=-1, sign_eq=0 means sign_j=+1, so term = -J.
        // Chunk0: sum of -1,-2,-3,-4 = -10
        // Chunk1: previous -10 + (-5-6) = -21
        // The XNOR inputs must be ~J because sign_xi=-1.
        // --------------------------------------------------------------------
        $display("Test 4: Two chunks, sign_xi=-1 (all sign_eq=0)");
        xnor_c0_neg = 0;
        xnor_c1_neg = 0;
        xnor_c0_neg |= pack_xnor_J(0, ~4'sd1) | pack_xnor_J(1, ~4'sd2) |
                       pack_xnor_J(2, ~4'sd3) | pack_xnor_J(3, ~4'sd4);
        xnor_c1_neg |= pack_xnor_J(0, ~4'sd5) | pack_xnor_J(1, ~4'sd6);
        run_test("Test4a chunk0", 0, 0, 4, 0, xnor_c0_neg, 4'b0000, -10);
        run_test("Test4b chunk1", 1, 1, 2, 0, xnor_c1_neg, 2'b00, -21);

        // --------------------------------------------------------------------
        // Test 5: Partial chunk, sign_xi=+1, J=[7,-2], valid=2, sign_eq=11
        // Sum = 7 + (-2) = 5
        // --------------------------------------------------------------------
        $display("Test 5: Partial chunk, sign_xi=+1");
        xnor_partial = 0;
        xnor_partial |= pack_xnor_J(0, 4'sd7);
        xnor_partial |= pack_xnor_J(1, -4'sd2);
        run_test("Test5", 0, 1, 2, 1, xnor_partial, 2'b11, 5);

        // --------------------------------------------------------------------
        // Summary
        // --------------------------------------------------------------------
        if (error_count == 0)
            $display("=== All tests PASSED ===");
        else begin
            $error("=== %0d tests FAILED ===", error_count);
            $finish(1);
        end
        $finish;
    end

    initial begin
        $dumpfile("tb_dotprod_phase2_sparse.vcd");
        $dumpvars(0, tb_dotprod_phase2_sparse);
    end

endmodule
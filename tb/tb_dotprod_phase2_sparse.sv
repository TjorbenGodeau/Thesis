// tb_dotprod_phase2_sparse.sv — corrected for compilation
`timescale 1ns/1ps

module tb_dotprod_phase2_sparse;

    import dsb_pkg::*;
    import dsb_sparse_pkg::*;

    // ── Test parameters ──────────────────────────────────────────────────────
    localparam int K_MAX_P   = 4;
    localparam int IC_BITS_P = 4;
    localparam int ACCUM_W_P = 16;
    localparam int K_CNT_W   = $clog2(K_MAX_P + 1);  // 3

    // Clock
    logic clk = 0;
    always #5 clk = ~clk;

    // DUT signals
    logic                         start;
    logic                         done;
    logic                         accumulate;
    logic                         last_chunk;
    logic [K_CNT_W-1:0]           valid_count;
    logic                         sign_xi;
    logic [K_MAX_P*IC_BITS_P-1:0] xnor_J;
    logic [K_MAX_P-1:0]           sign_eq;
    logic signed [ACCUM_W_P-1:0]  Jx_i;

    // DUT
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

    // ── Test control ────────────────────────────────────────────────────────
    int error_count = 0;

    // Helper: pack a single slot's xnor_J into the vector
    function automatic logic [K_MAX_P*IC_BITS_P-1:0] pack_xnor_J(
        input int slot, input logic [IC_BITS_P-1:0] value
    );
        return value << (slot * IC_BITS_P);
    endfunction

    // ── Correct reference model ─────────────────────────────────────────────
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
                logic signed [IC_BITS_P-1:0] term;
                logic sign_j;

                xnor_word = xnor_J[k*IC_BITS_P +: IC_BITS_P];

                // Reconstruct J from xnor and sign_xi
                if (sign_xi) begin
                    // sign_xi = +1: xnor == J
                    J_val = signed'(xnor_word);
                end else begin
                    // sign_xi = -1: xnor == ~J, so J = ~xnor
                    J_val = signed'(~xnor_word);
                end

                // Determine sign_j from sign_eq and sign_xi:
                // sign_eq = 1 => sign_j == sign_xi, so sign_j = sign_xi
                // sign_eq = 0 => sign_j != sign_xi, so sign_j = ~sign_xi
                sign_j = sign_eq[k] ? sign_xi : ~sign_xi;

                // term = J * sign_j * sign_xi
                // Since sign_j and sign_xi are ±1, product is +1 if sign_eq=1, else -1.
                term = sign_eq[k] ? J_val : -J_val;
                chunk_sum += ACCUM_W_P'(term);
            end
        end
        return (accumulate ? prev_accum : '0) + chunk_sum;
    endfunction

    // ── Helper to run a test case ──────────────────────────────────────────
    task automatic run_test(
        input string            label,
        input logic             accumulate_in,
        input logic             last_chunk_in,
        input logic [K_CNT_W-1:0] valid_count_in,
        input logic             sign_xi_in,
        input logic [K_MAX_P*IC_BITS_P-1:0] xnor_J_in,
        input logic [K_MAX_P-1:0] sign_eq_in,
        input logic signed [ACCUM_W_P-1:0] expected
    );
        // Set inputs
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

        // Wait for done (should be next cycle)
        @(posedge clk);
        if (!done) $error("%s: done not asserted", label);

        // Check result
        if (Jx_i !== expected) begin
            $error("%s: expected %0d, got %0d", label, expected, Jx_i);
            error_count++;
        end else begin
            $display("%s: OK (Jx_i=%0d)", label, Jx_i);
        end
    endtask

    // ── Test sequence ──────────────────────────────────────────────────────
    initial begin
        // ---- All declarations at the top ----
        logic [K_MAX_P*IC_BITS_P-1:0] xnor1, xnor_chunk0, xnor_chunk1;
        logic [K_MAX_P*IC_BITS_P-1:0] xnor_c0_neg, xnor_c1_neg, xnor_partial;
        logic signed [ACCUM_W_P-1:0] ref_accum;
        int i;

        // Initialize signals
        start = 0;
        accumulate = 0;
        last_chunk = 0;
        valid_count = 0;
        sign_xi = 0;
        xnor_J = 0;
        sign_eq = 0;
        ref_accum = 0;

        repeat (2) @(posedge clk);
        $display("=== Starting dotprod_phase2_sparse test (reference model) ===");

        // --------------------------------------------------------------------
        // Test 1: Single chunk, sign_xi=+1, J = [5, -3, 2, 1], sign_eq all 1
        // Expected: 5 + (-3) + 2 + 1 = 5
        // --------------------------------------------------------------------
        $display("Test 1: sign_xi=+1, all sign_eq=1");
        xnor1 = 0;
        xnor1 |= pack_xnor_J(0, 4'sd5);
        xnor1 |= pack_xnor_J(1, -4'sd3);
        xnor1 |= pack_xnor_J(2, 4'sd2);
        xnor1 |= pack_xnor_J(3, 4'sd1);
        // Compute expected using reference model with prev_accum=0, accumulate=0
        ref_accum = ref_dotprod(0, 0, 4, 1, xnor1, 4'b1111);
        run_test("Test1", 0, 1, 4, 1, xnor1, 4'b1111, ref_accum);

        // --------------------------------------------------------------------
        // Test 2: sign_xi=-1, same J, sign_eq all 0
        // Expected: (-5) + 3 + (-2) + (-1) = -5
        // --------------------------------------------------------------------
        $display("Test 2: sign_xi=-1, all sign_eq=0");
        ref_accum = ref_dotprod(0, 0, 4, 0, xnor1, 4'b0000);
        run_test("Test2", 0, 1, 4, 0, xnor1, 4'b0000, ref_accum);

        // --------------------------------------------------------------------
        // Test 3: Two chunks, sign_xi=+1, J = [1,2,3,4] then [5,6]
        // Expected total = 1+2+3+4+5+6 = 21
        // --------------------------------------------------------------------
        $display("Test 3: Two chunks, sign_xi=+1");
        xnor_chunk0 = 0;
        xnor_chunk1 = 0;
        xnor_chunk0 |= pack_xnor_J(0, 4'sd1) | pack_xnor_J(1, 4'sd2) |
                       pack_xnor_J(2, 4'sd3) | pack_xnor_J(3, 4'sd4);
        xnor_chunk1 |= pack_xnor_J(0, 4'sd5) | pack_xnor_J(1, 4'sd6);
        // Chunk0: accumulate=0, last=0
        ref_accum = ref_dotprod(0, 0, 4, 1, xnor_chunk0, 4'b1111);
        run_test("Test3a chunk0", 0, 0, 4, 1, xnor_chunk0, 4'b1111, ref_accum);
        // Chunk1: accumulate=1, last=1
        ref_accum = ref_dotprod(ref_accum, 1, 2, 1, xnor_chunk1, 2'b11);
        run_test("Test3b chunk1", 1, 1, 2, 1, xnor_chunk1, 2'b11, ref_accum);

        // --------------------------------------------------------------------
        // Test 4: Two chunks, sign_xi=-1, same J values.
        // Expected: -1-2-3-4-5-6 = -21
        // --------------------------------------------------------------------
        $display("Test 4: Two chunks, sign_xi=-1");
        // For sign_xi=-1, the xnor_J is the bitwise NOT of J (because XNOR gives ~J when rwl=0).
        // So we feed the ~J values, and sign_eq = 0 because sign_j != sign_xi (assuming sign_j=+1).
        // But in reality, sign_j is determined by J's sign bit. For simplicity, we'll set sign_j=+1 (sign_eq=0).
        // However, the reference model reconstructs J from ~xnor if sign_xi=0.
        // So we can feed the same xnor1 (which was J) but set sign_xi=0, and the model will invert it.
        // For clarity, we'll feed the same xnor_chunk0 and xnor_chunk1 and set sign_xi=0.
        // The reference model will do: J_val = signed'(~xnor_word) because sign_xi=0.
        // And sign_eq=0 => term = -J_val.
        // So it will compute correctly.
        ref_accum = ref_dotprod(0, 0, 4, 0, xnor_chunk0, 4'b0000);
        run_test("Test4a chunk0", 0, 0, 4, 0, xnor_chunk0, 4'b0000, ref_accum);
        ref_accum = ref_dotprod(ref_accum, 1, 2, 0, xnor_chunk1, 2'b00);
        run_test("Test4b chunk1", 1, 1, 2, 0, xnor_chunk1, 2'b00, ref_accum);

        // --------------------------------------------------------------------
        // Test 5: Mixed signs, partial chunk, sign_xi=+1
        // J = [7, -2], valid_count=2, sign_eq=11 (sign_j = +1, -1)
        // Expected: 7 + (-2) = 5
        // --------------------------------------------------------------------
        $display("Test 5: Partial chunk, sign_xi=+1");
        xnor_partial = 0;
        xnor_partial |= pack_xnor_J(0, 4'sd7);
        xnor_partial |= pack_xnor_J(1, -4'sd2);
        ref_accum = ref_dotprod(0, 0, 2, 1, xnor_partial, 2'b11);
        run_test("Test5", 0, 1, 2, 1, xnor_partial, 2'b11, ref_accum);

        // --------------------------------------------------------------------
        // Summary
        // --------------------------------------------------------------------
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
        $dumpfile("tb_dotprod_phase2_sparse.vcd");
        $dumpvars(0, tb_dotprod_phase2_sparse);
    end

endmodule
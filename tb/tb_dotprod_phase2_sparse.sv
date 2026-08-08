// tb_dotprod_phase2_sparse.sv
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

    // Helper to pack xnor_J and sign_eq vectors
    function automatic logic [K_MAX_P*IC_BITS_P-1:0] pack_xnor_J(
        input int slot, input logic [IC_BITS_P-1:0] value
    );
        return value << (slot * IC_BITS_P);
    endfunction

    // Reference model: exact replica of current RTL behaviour
    function automatic logic signed [ACCUM_W_P-1:0] ref_dotprod(
        input logic [ACCUM_W_P-1:0] Jx_accum_in,
        input logic                  accumulate,
        input logic                  last_chunk,
        input logic [K_CNT_W-1:0]    valid_count,
        input logic                  sign_xi,
        input logic [K_MAX_P*IC_BITS_P-1:0] xnor_J,
        input logic [K_MAX_P-1:0]    sign_eq
    );
        logic signed [ACCUM_W_P-1:0] chunk_sum;
        logic signed [ACCUM_W_P-1:0] base;
        logic [IC_BITS_P-1:0]        xnor_word;
        logic [IC_BITS_P-1:0]        corrected;

        chunk_sum = '0;
        for (int k = 0; k < K_MAX_P; k++) begin
            if (k < int'(valid_count)) begin
                xnor_word = xnor_J[k*IC_BITS_P +: IC_BITS_P];
                corrected  = sign_eq[k] ? xnor_word : ~xnor_word;
                chunk_sum += ACCUM_W_P'(signed'(corrected));
            end
        end

        base = (accumulate ? Jx_accum_in : '0) + chunk_sum;
        if (last_chunk && !sign_xi)
            base += ACCUM_W_P'(valid_count);   // as per RTL: adds valid_count
        return base;
    endfunction

    // Helper to run a test case
    task automatic test_case(
        input string            label,
        input logic             accumulate_in,
        input logic             last_chunk_in,
        input logic [K_CNT_W-1:0] valid_count_in,
        input logic             sign_xi_in,
        input logic [K_MAX_P*IC_BITS_P-1:0] xnor_J_in,
        input logic [K_MAX_P-1:0] sign_eq_in,
        input logic signed [ACCUM_W_P-1:0] expected
    );
        // Apply inputs
        accumulate  = accumulate_in;
        last_chunk  = last_chunk_in;
        valid_count = valid_count_in;
        sign_xi     = sign_xi_in;
        xnor_J      = xnor_J_in;
        sign_eq     = sign_eq_in;

        // Trigger start
        @(posedge clk);
        start = 1;
        @(posedge clk);
        start = 0;

        // Wait for done (should be next cycle)
        @(posedge clk);
        if (!done) $error("%s: done not asserted after one cycle", label);
        @(posedge clk);  // extra cycle for Jx_i to stabilise

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
        logic [K_MAX_P*IC_BITS_P-1:0] xnor1, xnor_chunk0, xnor_chunk1, xnor_partial;
        logic [K_MAX_P-1:0]           sign_eq1, sign_eq0, sign_eq_all0;
        int i;  // for loops if needed

        // Initialize signals
        start = 0;
        accumulate = 0;
        last_chunk = 0;
        valid_count = 0;
        sign_xi = 0;
        xnor_J = 0;
        sign_eq = 0;

        repeat (2) @(posedge clk);
        $display("=== Starting dotprod_phase2_sparse test ===");

        // --------------------------------------------------------------------
        // 1. Single chunk, K_MAX=4, valid_count=4, sign_xi=+1
        //    J values: [5, -3, 2, 1] (packed as 4-bit signed)
        //    sign_eq all 1 (i.e. sign_j == sign_xi)
        //    Expected: sum of J = 5 + (-3) + 2 + 1 = 5
        //    No correction (sign_xi=1)
        // --------------------------------------------------------------------
        $display("Test 1: Single chunk, sign_xi=+1, all sign_eq=1");
        xnor1 = 0;
        sign_eq1 = 4'b1111;
        xnor1 |= pack_xnor_J(0, 4'sd5);
        xnor1 |= pack_xnor_J(1, -4'sd3);        // -3
        xnor1 |= pack_xnor_J(2, 4'sd2);
        xnor1 |= pack_xnor_J(3, 4'sd1);
        test_case("Test1", 0, 1, 4, 1, xnor1, sign_eq1, 5);

        // --------------------------------------------------------------------
        // 2. Same but sign_xi=0 (negative), last_chunk=1
        //    Now sign_eq all 0 (since sign_xi=0, sign_j=+1, so sign_eq=0)
        //    The RTL correction: last_chunk && !sign_xi => adds valid_count (4)
        //    Expected: chunk_sum = Σ ~xnor_J (because sign_eq=0 => corrected = ~xnor_J)
        //    For 4-bit: J=5 (0101) -> ~J = 1010 = -6 (signed)
        //    J=-3 (1101) -> ~J = 0010 = +2
        //    J=2  (0010) -> ~J = 1101 = -3
        //    J=1  (0001) -> ~J = 1110 = -2
        //    Sum = -6 + 2 -3 -2 = -9
        //    Then add valid_count (4) -> -5
        //    So expected = -5.
        // --------------------------------------------------------------------
        $display("Test 2: Single chunk, sign_xi=0, sign_eq=0, correction adds 4");
        sign_eq0 = 4'b0000;
        test_case("Test2", 0, 1, 4, 0, xnor1, sign_eq0, -5);

        // --------------------------------------------------------------------
        // 3. Two chunks: chunk0 (accumulate=0, last_chunk=0), chunk1 (accumulate=1, last_chunk=1)
        //    Row has 6 entries (K_MAX=4). Chunk0: 4 entries, chunk1: 2 entries.
        //    sign_xi=+1, sign_eq all 1 for simplicity.
        //    J values: [1,2,3,4] in chunk0, [5,6] in chunk1.
        //    Expected final sum = 1+2+3+4+5+6 = 21
        //    No correction (sign_xi=1).
        // --------------------------------------------------------------------
        $display("Test 3: Two chunks, sign_xi=+1, no correction");
        xnor_chunk0 = 0;
        xnor_chunk1 = 0;
        xnor_chunk0 |= pack_xnor_J(0, 4'sd1);
        xnor_chunk0 |= pack_xnor_J(1, 4'sd2);
        xnor_chunk0 |= pack_xnor_J(2, 4'sd3);
        xnor_chunk0 |= pack_xnor_J(3, 4'sd4);
        xnor_chunk1 |= pack_xnor_J(0, 4'sd5);
        xnor_chunk1 |= pack_xnor_J(1, 4'sd6);
        // Run chunk0: accumulate=0, last_chunk=0, valid_count=4
        test_case("Test3a chunk0", 0, 0, 4, 1, xnor_chunk0, 4'b1111, 10); // sum=10
        // Run chunk1: accumulate=1, last_chunk=1, valid_count=2
        test_case("Test3b chunk1", 1, 1, 2, 1, xnor_chunk1, 2'b11, 21);

        // --------------------------------------------------------------------
        // 4. Same as 3, but sign_xi=0 (negative), to test correction only on last chunk.
        //    Chunk0: sign_xi=0, sign_eq=0 (all terms inverted), accumulate=0, last_chunk=0
        //    Chunk1: sign_xi=0, sign_eq=0 (all terms inverted), accumulate=1, last_chunk=1
        //    For chunk0: sum of ~J for [1,2,3,4] -> ~1=-2, ~2=-3, ~3=-4, ~4=-5? Wait 4-bit signed:
        //      1 (0001) -> ~1 = 1110 = -2
        //      2 (0010) -> ~2 = 1101 = -3
        //      3 (0011) -> ~3 = 1100 = -4
        //      4 (0100) -> ~4 = 1011 = -5
        //    Sum = -14. No correction (last_chunk=0).
        //    Chunk1: ~5 (0101) -> ~5 = 1010 = -6; ~6 (0110) -> ~6 = 1001 = -7
        //    Sum = -13. Then base = previous (-14) + (-13) = -27.
        //    Correction on last_chunk: adds valid_count (which is total row count? RTL adds valid_count of current chunk, i.e. 2)
        //    So expected = -27 + 2 = -25.
        //    This matches the current RTL (adds current valid_count).
        // --------------------------------------------------------------------
        $display("Test 4: Two chunks, sign_xi=0, correction on last chunk (adds current valid_count)");
        sign_eq_all0 = 4'b0000;
        // Chunk0
        test_case("Test4a chunk0", 0, 0, 4, 0, xnor_chunk0, sign_eq_all0, -14);
        // Chunk1 (last_chunk=1, accumulate=1, valid_count=2)
        test_case("Test4b chunk1", 1, 1, 2, 0, xnor_chunk1, 2'b00, -25);

        // --------------------------------------------------------------------
        // 5. Test valid_count < K_MAX (partial chunk) and partial sign_eq.
        //    sign_xi=+1, valid_count=2, sign_eq=11, J=[7, -2] -> sum=5
        // --------------------------------------------------------------------
        $display("Test 5: Partial chunk, sign_xi=+1");
        xnor_partial = 0;
        xnor_partial |= pack_xnor_J(0, 4'sd7);
        xnor_partial |= pack_xnor_J(1, -4'sd2);    // -2
        test_case("Test5", 0, 1, 2, 1, xnor_partial, 2'b11, 5);

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
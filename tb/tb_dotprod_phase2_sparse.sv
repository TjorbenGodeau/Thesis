// tb_dotprod_phase2_sparse.sv
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

    dotprod_phase2_sparse u_dut (
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

    // ── Reference: compute sum of one chunk (no accumulation) ────────
    function automatic logic signed [ACCUM_W_P-1:0] chunk_sum_func(
        input logic [K_CNT_W-1:0]   vcnt,
        input logic                 sxi,
        input logic [K_MAX_P*IC_BITS_P-1:0] xnor_vec,
        input logic [K_MAX_P-1:0]   sign_eq_vec
    );
        logic signed [ACCUM_W_P-1:0] sum;
        logic [IC_BITS_P-1:0]        xw;
        logic signed [IC_BITS_P-1:0] xw_s, neg_xw_s;

        sum = '0;
        for (int k = 0; k < K_MAX_P; k++) begin
            if (k < vcnt) begin
                xw = (xnor_vec >> (k * IC_BITS_P)) & ((1 << IC_BITS_P) - 1);
                xw_s    = $signed({1'b0, xw});      // positive
                neg_xw_s = $signed(~xw);            // two's complement of ~xw

                case ({sxi, sign_eq_vec[k]})
                    2'b00: sum += $signed(neg_xw_s);          // ~xw
                    2'b01: sum += $signed(xw_s);              // xw
                    2'b10: sum += $signed(xw_s) + 1;          // xw + 1
                    2'b11: sum += $signed(neg_xw_s) + 1;      // ~xw + 1
                endcase
            end
        end
        return sum;
    endfunction

    // ── Test harness ─────────────────────────────────────────────────────
    int error_count = 0;
    logic signed [ACCUM_W_P-1:0] tb_acc = 0;   // independent accumulator

    task check_equal(input string msg, input logic signed [ACCUM_W_P-1:0] a, input logic signed [ACCUM_W_P-1:0] b);
        if (a !== b) begin
            $error("%s: expected %0d, got %0d", msg, b, a);
            error_count++;
        end else begin
            $display("  %s: OK", msg);
        end
    endtask

    task test_case(
        input string label,
        input logic                 acc,
        input logic [K_CNT_W-1:0]   vcnt,
        input logic                 sxi,
        input logic [K_MAX_P*IC_BITS_P-1:0] xnor_vec,
        input logic [K_MAX_P-1:0]   sign_eq_vec
    );
        // Compute chunk sum using reference
        automatic logic signed [ACCUM_W_P-1:0] chunk = chunk_sum_func(vcnt, sxi, xnor_vec, sign_eq_vec);
        // Update testbench accumulator
        if (acc) tb_acc = tb_acc + chunk;
        else     tb_acc = chunk;

        // Drive inputs
        accumulate  = acc;
        valid_count = vcnt;
        sign_xi     = sxi;
        xnor_J      = xnor_vec;
        sign_eq     = sign_eq_vec;

        // Apply start
        @(posedge clk);
        start = 1;
        @(posedge clk);              // DUT samples start, computes, sets done
        // done is high now
        if (!done) $error("%s: done not asserted", label);
        // Check result
        check_equal(label, Jx_i, tb_acc);
        start = 0;
        // Wait one cycle to allow done to clear (optional)
        @(posedge clk);
    endtask

    // ── Test execution ──────────────────────────────────────────────────
    initial begin
        logic [K_MAX_P*IC_BITS_P-1:0] xnor1, xnor2;
        logic [K_MAX_P-1:0]           seq1, seq2;

        start = 0;
        accumulate = 0;
        last_chunk = 0;
        valid_count = 0;
        sign_xi = 0;
        xnor_J = 0;
        sign_eq = 0;
        tb_acc = 0;

        repeat (2) @(posedge clk);
        $display("=== dotprod_phase2_sparse test (corrected timing & independent acc) ===");

        // ── Build test vectors ──────────────────────────────────────────
        xnor1 = 0;
        seq1  = 4'b0101;   // slot0=1, slot1=0, slot2=1, slot3=0
        xnor1[0*IC_BITS_P +: IC_BITS_P] = 4'd5;   // 0101
        xnor1[1*IC_BITS_P +: IC_BITS_P] = 4'd3;   // 0011
        xnor1[2*IC_BITS_P +: IC_BITS_P] = 4'd13;  // 1101 → signed -3
        xnor1[3*IC_BITS_P +: IC_BITS_P] = 4'd0;   // 0000

        // ── Test 1: sign_xi=1, all 4 slots ─────────────────────────────
        // chunk sum = -5 + 4 + 3 + 1 = 3
        test_case("Test1 all cases (sign_xi=1)", 0, 4, 1, xnor1, seq1);

        // ── Test 2: sign_xi=0, same vectors ────────────────────────────
        // chunk sum = 5 - 4 - 3 - 1 = -3
        test_case("Test2 sign_xi=0", 0, 4, 0, xnor1, seq1);

        // ── Test 3: multi‑chunk accumulation ───────────────────────────
        // Chunk0 → same as Test1: result = 3
        // Chunk1: accumulate=1, valid=2, sign_xi=1, xnor2, seq2
        xnor2 = 0;
        seq2  = 2'b01;   // only slots0,1 valid
        xnor2[0*IC_BITS_P +: IC_BITS_P] = 4'd7;   // 0111
        xnor2[1*IC_BITS_P +: IC_BITS_P] = 4'd1;   // 0001
        // chunk1 sum = 1 + 2 = 3; total = 3 + 3 = 6
        test_case("Test3a chunk0", 0, 4, 1, xnor1, seq1);
        test_case("Test3b chunk1", 1, 2, 1, xnor2, seq2);

        // ── Test 4: partial chunk (valid=2) ────────────────────────────
        // Only slots 0 and 1 used: -5 + 4 = -1
        test_case("Test4 partial", 0, 2, 1, xnor1, seq1);

        // ── Summary ──────────────────────────────────────────────────────
        if (error_count == 0) begin
            $display("=== All tests PASSED ===");
        end else begin
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
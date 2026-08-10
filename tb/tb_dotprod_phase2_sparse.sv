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

    // ── Reference model (exact copy of RTL case table) ──────────────────
    function automatic logic signed [ACCUM_W_P-1:0] ref_dotprod(
        input logic [ACCUM_W_P-1:0] prev,
        input logic                 acc,
        input logic [K_CNT_W-1:0]   vcnt,
        input logic                 sxi,
        input logic [K_MAX_P*IC_BITS_P-1:0] xnor,
        input logic [K_MAX_P-1:0]   seq
    );
        logic signed [ACCUM_W_P-1:0] sum = '0;
        for (int k=0; k<K_MAX_P; k++) begin
            if (k < int'(vcnt)) begin
                logic [IC_BITS_P-1:0] xw = xnor[k*IC_BITS_P +: IC_BITS_P];
                logic signed [ACCUM_W_P-1:0] term;
                case ({sxi, seq[k]})
                    2'b00: term = ACCUM_W_P'(signed'(~xw));
                    2'b01: term = ACCUM_W_P'(signed'(xw) + 1);
                    2'b10: term = ACCUM_W_P'(signed'(~xw) + 1);
                    2'b11: term = ACCUM_W_P'(signed'(xw));
                endcase
                sum += term;
            end
        end
        return (acc ? prev : '0) + sum;
    endfunction

    // ── Test harness ─────────────────────────────────────────────────────
    int error_count = 0;

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
        input logic [K_MAX_P*IC_BITS_P-1:0] xnor,
        input logic [K_MAX_P-1:0]   seq
    );
        automatic logic signed [ACCUM_W_P-1:0] expected;
        expected = ref_dotprod(Jx_i, acc, vcnt, sxi, xnor, seq);
        accumulate  = acc;
        valid_count = vcnt;
        sign_xi     = sxi;
        xnor_J      = xnor;
        sign_eq     = seq;
        @(posedge clk);
        start = 1;
        @(posedge clk);
        start = 0;
        @(posedge clk);
        if (!done) $error("%s: done not asserted", label);
        @(posedge clk);
        check_equal(label, Jx_i, expected);
    endtask

    initial begin
        start = 0;
        accumulate = 0;
        last_chunk = 0;
        valid_count = 0;
        sign_xi = 0;
        xnor_J = 0;
        sign_eq = 0;

        repeat (2) @(posedge clk);
        $display("=== dotprod_phase2_sparse test (case table) ===");

        // ── 1. Single chunk, all four combinations ──────────────────────
        // We set xnor_J and sign_eq to exercise each case.
        // For each slot we choose a base value and compute expected manually.
        // Slot0: sign_xi=1, sign_eq=1 -> term = xnor
        // Slot1: sign_xi=1, sign_eq=0 -> term = ~xnor + 1
        // Slot2: sign_xi=0, sign_eq=1 -> term = xnor + 1
        // Slot3: sign_xi=0, sign_eq=0 -> term = ~xnor
        logic [K_MAX_P*IC_BITS_P-1:0] xnor1 = 0;
        logic [K_MAX_P-1:0]           seq1   = 4'b0101; // slot0=1, slot1=0, slot2=1, slot3=0
        xnor1[0*IC_BITS_P +: IC_BITS_P] = 4'd5;   // case 11 -> term=5
        xnor1[1*IC_BITS_P +: IC_BITS_P] = 4'd3;   // case 10 -> ~3+1 = -4+1 = -3
        xnor1[2*IC_BITS_P +: IC_BITS_P] = 4'd-3;  // case 01 -> -3+1 = -2
        xnor1[3*IC_BITS_P +: IC_BITS_P] = 4'd0;   // case 00 -> ~0 = -1
        // Sum = 5 + (-3) + (-2) + (-1) = -1
        test_case("Test1 all cases (sign_xi=1)", 0, 4, 1, xnor1, seq1);

        // ── 2. Sign_xi=0 with same xnor/seq ────────────────────────────
        // Cases change:
        // slot0: sign_xi=0, seq=1 -> case 01: xnor+1 = 5+1=6
        // slot1: sign_xi=0, seq=0 -> case 00: ~xnor = ~3 = -4
        // slot2: sign_xi=0, seq=1 -> case 01: -3+1 = -2
        // slot3: sign_xi=0, seq=0 -> case 00: ~0 = -1
        // Sum = 6 + (-4) + (-2) + (-1) = -1
        test_case("Test2 sign_xi=0", 0, 4, 0, xnor1, seq1);

        // ── 3. Multi‑chunk accumulation ─────────────────────────────────
        // Chunk0: accumulate=0, result = -1
        // Chunk1: accumulate=1, valid=2, sign_xi=1, xnor2, seq2
        logic [K_MAX_P*IC_BITS_P-1:0] xnor2 = 0;
        logic [K_MAX_P-1:0]           seq2   = 2'b01; // only slots0,1 valid
        xnor2[0*IC_BITS_P +: IC_BITS_P] = 4'd7;   // case 11 -> term=7
        xnor2[1*IC_BITS_P +: IC_BITS_P] = 4'd1;   // case 10 -> ~1+1 = -2+1 = -1
        // Sum chunk1 = 7 + (-1) = 6, total = -1 + 6 = 5
        test_case("Test3a chunk0", 0, 4, 1, xnor1, seq1);
        test_case("Test3b chunk1", 1, 2, 1, xnor2, seq2);

        // ── 4. Partial chunk (valid_count < K_MAX) ──────────────────────
        // Only slots 0 and 1 used: sum = 5 + (-3) = 2
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
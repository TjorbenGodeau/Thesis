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
    logic [K_MAX_P*IC_BITS_P-1:0] xnor_J;          // DUT port name remains
    logic [K_MAX_P-1:0]           sign_eq;         // DUT port name remains
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

    // ── Reference model (matches RTL case table) ──────────────────────
    function automatic logic signed [ACCUM_W_P-1:0] ref_dotprod(
        input logic [ACCUM_W_P-1:0] prev,
        input logic                 acc,
        input logic [K_CNT_W-1:0]   vcnt,
        input logic                 sxi,
        input logic [K_MAX_P*IC_BITS_P-1:0] xnor_vec,   // renamed to avoid keyword
        input logic [K_MAX_P-1:0]   sign_eq_vec        // renamed for clarity
    );
        logic signed [ACCUM_W_P-1:0] sum;
        logic [IC_BITS_P-1:0]        xw;
        logic signed [IC_BITS_P-1:0] xw_s, neg_xw_s;

        sum = '0;
        for (int k = 0; k < K_MAX_P; k++) begin
            if (k < vcnt) begin
                // Extract the k‑th field using shift & mask
                xw = (xnor_vec >> (k * IC_BITS_P)) & ((1 << IC_BITS_P) - 1);
                xw_s    = $signed({1'b0, xw});        // convert to signed (positive)
                neg_xw_s = $signed(~xw);              // two's complement of bitwise NOT

                case ({sxi, sign_eq_vec[k]})
                    2'b00: sum += $signed(neg_xw_s);          // ~xw
                    2'b01: sum += $signed(xw_s);              // xw
                    2'b10: sum += $signed(xw_s) + 1;          // xw + 1
                    2'b11: sum += $signed(neg_xw_s) + 1;      // ~xw + 1
                endcase
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
        input logic [K_MAX_P*IC_BITS_P-1:0] xnor_vec,
        input logic [K_MAX_P-1:0]   sign_eq_vec
    );
        automatic logic signed [ACCUM_W_P-1:0] expected;
        expected = ref_dotprod(Jx_i, acc, vcnt, sxi, xnor_vec, sign_eq_vec);
        accumulate  = acc;
        valid_count = vcnt;
        sign_xi     = sxi;
        xnor_J      = xnor_vec;
        sign_eq     = sign_eq_vec;
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
        $display("=== dotprod_phase2_sparse test (corrected) ===");

        // ── Test vectors ────────────────────────────────────────────────
        logic [K_MAX_P*IC_BITS_P-1:0] xnor1 = 0;
        logic [K_MAX_P-1:0]           seq1  = 4'b0101;   // slot0=1, slot1=0, slot2=1, slot3=0

        // slot0: xw = 5  (4'b0101)
        xnor1[0*IC_BITS_P +: IC_BITS_P] = 4'd5;
        // slot1: xw = 3  (4'b0011)
        xnor1[1*IC_BITS_P +: IC_BITS_P] = 4'd3;
        // slot2: xw = 13 (4'b1101) → signed -3
        xnor1[2*IC_BITS_P +: IC_BITS_P] = 4'd13;
        // slot3: xw = 0
        xnor1[3*IC_BITS_P +: IC_BITS_P] = 4'd0;

        // ── Test 1: sign_xi=1, all 4 slots ─────────────────────────────
        // Expected sum:
        // slot0: sxi=1, seq=1 → case 2'b11 → ~5+1 = -6+1 = -5
        // slot1: sxi=1, seq=0 → case 2'b10 → 3+1   =  4
        // slot2: sxi=1, seq=1 → case 2'b11 → ~13+1 =  2+1 =  3   (since ~13=2)
        // slot3: sxi=1, seq=0 → case 2'b10 → 0+1    =  1
        // Sum = -5 + 4 + 3 + 1 = 3
        test_case("Test1 all cases (sign_xi=1)", 0, 4, 1, xnor1, seq1);

        // ── Test 2: sign_xi=0, same vectors ────────────────────────────
        // slot0: sxi=0, seq=1 → case 2'b01 → xw = 5
        // slot1: sxi=0, seq=0 → case 2'b00 → ~xw = -4
        // slot2: sxi=0, seq=1 → case 2'b01 → xw = -3
        // slot3: sxi=0, seq=0 → case 2'b00 → ~xw = -1
        // Sum = 5 - 4 - 3 - 1 = -3
        test_case("Test2 sign_xi=0", 0, 4, 0, xnor1, seq1);

        // ── Test 3: multi‑chunk accumulation ───────────────────────────
        // Chunk0 → same as Test1: result = 3
        // Chunk1: accumulate=1, valid=2, sign_xi=1, xnor2, seq2
        logic [K_MAX_P*IC_BITS_P-1:0] xnor2 = 0;
        logic [K_MAX_P-1:0]           seq2  = 2'b01;   // only slots0,1 valid
        xnor2[0*IC_BITS_P +: IC_BITS_P] = 4'd7;   // xw=7
        xnor2[1*IC_BITS_P +: IC_BITS_P] = 4'd1;   // xw=1
        // slot0: sxi=1, seq=1 → case 2'b11 → ~7+1 = 0+1 = 1
        // slot1: sxi=1, seq=0 → case 2'b10 → 1+1 = 2
        // chunk sum = 1+2 = 3
        // total = 3 + 3 = 6
        test_case("Test3a chunk0", 0, 4, 1, xnor1, seq1);
        test_case("Test3b chunk1", 1, 2, 1, xnor2, seq2);

        // ── Test 4: partial chunk (valid=2) ────────────────────────────
        // Only slots 0 and 1 used:
        // slot0: -5, slot1: 4 → sum = -1
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
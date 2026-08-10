// tb_dotprod_phase2_sparse.sv
`timescale 1ns/1ps

module tb_dotprod_phase2_sparse;

    import dsb_pkg::*;
    import dsb_sparse_pkg::*;

    localparam int K_MAX_P   = 4;
    localparam int IC_BITS_P = 4;
    localparam int ACCUM_W_P = 16;
    localparam int K_CNT_W   = $clog2(K_MAX_P + 1);

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

    // Instantiate DUT
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

    // ===================================================================
    //  GOLDEN MODEL (per‑chunk sum, no accumulation)
    // ===================================================================
    function automatic logic signed [ACCUM_W_P-1:0] golden_chunk_sum(
        input logic [K_CNT_W-1:0]   vcnt,
        input logic                 sxi,
        input logic [K_MAX_P*IC_BITS_P-1:0] xnor_vec,
        input logic [K_MAX_P-1:0]   seq
    );
        logic signed [ACCUM_W_P-1:0] sum;
        logic [IC_BITS_P-1:0]        xw;
        logic signed [IC_BITS_P-1:0] xw_s, neg_xw_s;

        sum = '0;
        for (int k = 0; k < K_MAX_P; k++) begin
            if (k < vcnt) begin
                // Extract k‑th 4‑bit field
                xw = (xnor_vec >> (k*IC_BITS_P)) & ((1<<IC_BITS_P)-1);
                xw_s    = $signed({1'b0, xw});       // positive
                neg_xw_s = $signed(~xw);             // two's complement of ~xw

                // RTL case table: {sign_xi, sign_eq[k]}
                case ({sxi, seq[k]})
                    2'b00: sum += $signed(neg_xw_s);       // ~xw
                    2'b01: sum += $signed(xw_s);           //  xw
                    2'b10: sum += $signed(xw_s) + 1;       //  xw + 1
                    2'b11: sum += $signed(neg_xw_s) + 1;   // ~xw + 1
                endcase
            end
        end
        return sum;
    endfunction

    // ===================================================================
    //  TEST INFRASTRUCTURE
    // ===================================================================
    int error_count = 0;
    logic signed [ACCUM_W_P-1:0] expected_acc;   // our independent accumulator

    // Check and compare
    task check_result(input string label, input logic signed [ACCUM_W_P-1:0] actual);
        if (actual !== expected_acc) begin
            $error("%s: expected %0d, got %0d", label, expected_acc, actual);
            error_count++;
        end else begin
            $display("  %s: OK (result = %0d)", label, actual);
        end
    endtask

    // Single test case
    task run_test(
        input string label,
        input logic                 acc,          // accumulate flag for this chunk
        input logic [K_CNT_W-1:0]   vcnt,
        input logic                 sxi,
        input logic [K_MAX_P*IC_BITS_P-1:0] xnor_vec,
        input logic [K_MAX_P-1:0]   seq
    );
        // 1. Compute chunk sum using the golden model
        automatic logic signed [ACCUM_W_P-1:0] chunk = golden_chunk_sum(vcnt, sxi, xnor_vec, seq);

        // 2. Update the testbench accumulator (matches what DUT will do)
        if (acc) expected_acc += chunk;
        else     expected_acc  = chunk;

        // 3. Drive inputs
        accumulate  = acc;
        valid_count = vcnt;
        sign_xi     = sxi;
        xnor_J      = xnor_vec;
        sign_eq     = seq;

        // 4. Apply start pulse
        @(posedge clk);
        start = 1;
        @(posedge clk);               // DUT samples start and calculates
        start = 0;

        // 5. Check done and result
        if (!done) $error("%s: done not asserted", label);
        @(negedge clk);              // Give a little margin
        check_result(label, Jx_i);
        @(posedge clk);              // Let done clear (optional)
    endtask

    // ===================================================================
    //  TEST SEQUENCE
    // ===================================================================
    initial begin
        // Declare test vectors
        logic [K_MAX_P*IC_BITS_P-1:0] xnor1, xnor2;
        logic [K_MAX_P-1:0]           seq1, seq2;

        // Initialise
        start = 0;
        accumulate = 0;
        last_chunk = 0;
        valid_count = 0;
        sign_xi = 0;
        xnor_J = 0;
        sign_eq = 0;
        expected_acc = 0;

        repeat (2) @(posedge clk);
        $display("=== Clear, step‑by‑step test for dotprod_phase2_sparse ===");

        // ----------------------------------------------------------------
        // Build common test vector xnor1, seq1 (4 slots)
        // ----------------------------------------------------------------
        xnor1 = 0;
        seq1  = 4'b0101;   // slot0: seq=1, slot1: seq=0, slot2: seq=1, slot3: seq=0
        xnor1[0*IC_BITS_P +: IC_BITS_P] = 4'd5;   // 0101
        xnor1[1*IC_BITS_P +: IC_BITS_P] = 4'd3;   // 0011
        xnor1[2*IC_BITS_P +: IC_BITS_P] = 4'd13;  // 1101 (signed -3)
        xnor1[3*IC_BITS_P +: IC_BITS_P] = 4'd0;   // 0000

        // -----------------------------------------------------------------
        // Test 1: sign_xi=1, all 4 slots, accumulate=0 (start fresh)
        // Expected sum: 
        //   slot0: sxi=1, seq=1  → 2'b11 → ~5+1 = -6+1 = -5
        //   slot1: sxi=1, seq=0  → 2'b10 →  3+1 =  4
        //   slot2: sxi=1, seq=1  → 2'b11 → ~13+1 = 2+1 = 3  (since ~13=2)
        //   slot3: sxi=1, seq=0  → 2'b10 →  0+1 =  1
        // Total = -5 + 4 + 3 + 1 = 3
        // -----------------------------------------------------------------
        run_test("Test1: sign_xi=1, all slots", 0, 4, 1, xnor1, seq1);

        // -----------------------------------------------------------------
        // Test 2: sign_xi=0, same vectors, accumulate=0
        // Expected:
        //   slot0: sxi=0, seq=1 → 2'b01 →  5
        //   slot1: sxi=0, seq=0 → 2'b00 → ~3 = -4
        //   slot2: sxi=0, seq=1 → 2'b01 → -3
        //   slot3: sxi=0, seq=0 → 2'b00 → ~0 = -1
        // Total = 5 - 4 - 3 - 1 = -3
        // -----------------------------------------------------------------
        run_test("Test2: sign_xi=0, all slots", 0, 4, 0, xnor1, seq1);

        // -----------------------------------------------------------------
        // Test 3: Multi‑chunk accumulation
        // Chunk0: same as Test1 → result = 3 (accumulate=0)
        // Chunk1: accumulate=1, valid=2, sign_xi=1, new xnor2, seq2
        //   xnor2: slot0=7 (0111), slot1=1 (0001), only these two valid
        //   slot0: sxi=1, seq=1 → 2'b11 → ~7+1 = 0+1 = 1
        //   slot1: sxi=1, seq=0 → 2'b10 →  1+1 = 2
        //   chunk sum = 3
        //   total = 3 + 3 = 6
        // -----------------------------------------------------------------
        xnor2 = 0;
        seq2  = 2'b01;   // only slots 0,1 used
        xnor2[0*IC_BITS_P +: IC_BITS_P] = 4'd7;   // 0111
        xnor2[1*IC_BITS_P +: IC_BITS_P] = 4'd1;   // 0001

        run_test("Test3a: chunk0 (accum=0)", 0, 4, 1, xnor1, seq1);
        run_test("Test3b: chunk1 (accum=1)", 1, 2, 1, xnor2, seq2);

        // -----------------------------------------------------------------
        // Test 4: Partial chunk (valid=2) using xnor1, seq1
        // Only slots 0 and 1 used:
        //   slot0: -5, slot1: 4 → sum = -1
        // -----------------------------------------------------------------
        run_test("Test4: partial (valid=2)", 0, 2, 1, xnor1, seq1);

        // -----------------------------------------------------------------
        // Summary
        // -----------------------------------------------------------------
        if (error_count == 0) begin
            $display("=== ALL TESTS PASSED ===");
        end else begin
            $error("=== %0d TESTS FAILED ===", error_count);
            $finish(1);
        end
        $finish;
    end

    initial begin
        $dumpfile("tb_dotprod_phase2_sparse.vcd");
        $dumpvars(0, tb_dotprod_phase2_sparse);
    end

endmodule
// tb_sparse_compute_memory_advanced.sv
// Advanced self-checking testbench for sparse_compute_memory
`timescale 1ns/1ps

module tb_sparse_compute_memory_advanced;

    import dsb_pkg::*;
    import dsb_sparse_pkg::*;

    // Parameters (K_MAX=4, IC_BITS=4, COL_W=3)
    localparam int K_MAX_T  = 4;
    localparam int IC_BITS_T = 4;
    localparam int COL_W_T   = 3;

    // Clock
    logic clk = 0;
    always #5 clk = ~clk;

    // DUT signals
    logic                        wr_en;
    logic [$clog2(K_MAX_T)-1:0]  wr_slot;
    logic [COL_W_T-1:0]          wr_col_idx;
    logic [IC_BITS_T-1:0]        wr_J;
    logic                        wr_sign_j;

    logic                        rd_col_en;
    logic [$clog2(K_MAX_T)-1:0]  rd_slot;
    logic [COL_W_T-1:0]          rd_col_idx;
    logic                        rd_col_valid;

    logic                        precharge;
    logic [K_MAX_T-1:0]          rwl;
    logic [K_MAX_T*IC_BITS_T-1:0] rbl_J;
    logic [K_MAX_T-1:0]          rbl_signs;

    // DUT
    sparse_compute_memory #(
        .K_MAX_P  (K_MAX_T),
        .IC_BITS_P(IC_BITS_T),
        .COL_W_P  (COL_W_T)
    ) u_dut (
        .clk         (clk),
        .wr_en       (wr_en),
        .wr_slot     (wr_slot),
        .wr_col_idx  (wr_col_idx),
        .wr_J        (wr_J),
        .wr_sign_j   (wr_sign_j),
        .rd_col_en   (rd_col_en),
        .rd_slot     (rd_slot),
        .rd_col_idx  (rd_col_idx),
        .rd_col_valid(rd_col_valid),
        .precharge   (precharge),
        .rwl         (rwl),
        .rbl_J       (rbl_J),
        .rbl_signs   (rbl_signs)
    );

    // ─── Test control ────────────────────────────────────────────────────
    int error_count = 0;
    int test_pass_count = 0;

    // Helper: extract a slot's J bits from rbl_J
    function automatic logic [IC_BITS_T-1:0] get_rbl_J(input int slot);
        return rbl_J[slot*IC_BITS_T +: IC_BITS_T];
    endfunction

    // Helper: expected XNOR
    function automatic logic [IC_BITS_T-1:0] expected_xnor_J(input logic [IC_BITS_T-1:0] stored, input logic rwl);
        return (rwl ? stored : ~stored);
    endfunction
    function automatic logic expected_xnor_bit(input logic stored, input logic rwl);
        return ~(stored ^ rwl);
    endfunction

    // Check tasks
    task check_equal_vec(input string msg, input logic [IC_BITS_T-1:0] a, input logic [IC_BITS_T-1:0] b);
        if (a !== b) begin
            $error("%s: expected %h, got %h", msg, b, a);
            error_count++;
        end else begin
            $display("  %s: OK", msg);
            test_pass_count++;
        end
    endtask

    task check_equal_bit(input string msg, input logic a, input logic b);
        if (a !== b) begin
            $error("%s: expected %b, got %b", msg, b, a);
            error_count++;
        end else begin
            $display("  %s: OK", msg);
            test_pass_count++;
        end
    endtask

    // ─── Stored expected data ──────────────────────────────────────────
    logic [COL_W_T-1:0]          exp_col [K_MAX_T];
    logic [IC_BITS_T-1:0]        exp_J  [K_MAX_T];
    logic                        exp_sign[K_MAX_T];
    logic                        stored_sign_verified [K_MAX_T];

    // ─── Write task with verification ──────────────────────────────────
    task automatic write_slot_and_verify(   // ★ Make automatic
        input int slot,
        input logic [COL_W_T-1:0] col,
        input logic [IC_BITS_T-1:0] J,
        input logic sign
    );
        // ★ All declarations at the top of the task
        logic actual_sign;
        logic exp;
        logic [IC_BITS_T-1:0] actual_J;
        logic [IC_BITS_T-1:0] exp_J_after;

        // Write
        @(posedge clk);
        wr_en = 1; wr_slot = slot; wr_col_idx = col; wr_J = J; wr_sign_j = sign;
        exp_col[slot] = col; exp_J[slot] = J; exp_sign[slot] = sign;
        $display("[%0t] Writing slot%d: col=%d, J=%h, sign=%b", $time, slot, col, J, sign);
        @(posedge clk);
        wr_en = 0;

        // Immediately read back stored value using rwl=1 (XNOR with 1 gives stored)
        // Wait one cycle for combinational settle after changing rwl
        @(posedge clk);
        rwl = (1 << slot);   // only this slot gets rwl=1, others get 0
        #1;
        // Check sign
        actual_sign = rbl_signs[slot];
        exp = expected_xnor_bit(sign, 1); // rwl=1 => XNOR(sign,1) = sign
        check_equal_bit($sformatf("slot%d sign after write (rwl=1)", slot), actual_sign, exp);
        // Check J bits
        actual_J = get_rbl_J(slot);
        exp_J_after = expected_xnor_J(J, 1); // = J
        check_equal_vec($sformatf("slot%d J after write (rwl=1)", slot), actual_J, exp_J_after);

        // Also check other slots are unchanged (optional)
        // We'll do that in a separate test.

        // Also verify col_idx via read port
        @(posedge clk);
        rd_col_en = 1; rd_slot = slot;
        @(posedge clk);
        // Now rd_col_valid is high
        check_equal_vec($sformatf("col_idx slot%d", slot), rd_col_idx, col);
        if (rd_col_valid !== 1) $error("rd_col_valid not high for slot%d", slot);
        rd_col_en = 0;

        // Reset rwl to avoid affecting later tests
        @(posedge clk);
        rwl = 0;
    endtask

    // ─── Test sequence ──────────────────────────────────────────────────
    initial begin
        // Initialize
        wr_en = 0; rd_col_en = 0; precharge = 0; rwl = 0;
        for (int i=0; i<K_MAX_T; i++) begin
            exp_col[i]  = '0;
            exp_J[i]    = '0;
            exp_sign[i] = 1'b0;
        end

        repeat (2) @(posedge clk);
        $display("=== Starting advanced sparse_compute_memory test ===");

        // -----------------------------------------------------------------
        // 1. Write each slot and verify individually
        // -----------------------------------------------------------------
        $display("Test 1: Write each slot and verify stored values");
        write_slot_and_verify(0, 2, 4'b0101, 1'b1);
        write_slot_and_verify(1, 7, 4'b1010, 1'b0);
        write_slot_and_verify(2, 3, 4'b0111, 1'b1);
        write_slot_and_verify(3, 0, 4'b1111, 1'b0);

        // -----------------------------------------------------------------
        // 2. Compute with rwl == sign_j (all sign_eq should be 1)
        //    rwl = {slot3=0, slot2=1, slot1=0, slot0=1} => 4'b0101
        // -----------------------------------------------------------------
        $display("Test 2: Compute with rwl == sign_j");
        rwl = 4'b0101;   // slot0=1, slot1=0, slot2=1, slot3=0
        #1;
        for (int s=0; s<K_MAX_T; s++) begin
            // ★ All declarations first
            logic exp_s;
            logic [IC_BITS_T-1:0] exp_j;
            // Then assignments
            exp_s = expected_xnor_bit(exp_sign[s], rwl[s]);
            exp_j = expected_xnor_J(exp_J[s], rwl[s]);
            // Then checks
            check_equal_bit($sformatf("slot%d sign_eq", s), rbl_signs[s], exp_s);
            check_equal_vec($sformatf("slot%d J bits", s), get_rbl_J(s), exp_j);
        end

        // -----------------------------------------------------------------
        // 3. Compute with rwl != sign_j (all sign_eq should be 0)
        //    rwl = {slot3=1, slot2=0, slot1=1, slot0=0} => 4'b1010
        // -----------------------------------------------------------------
        $display("Test 3: Compute with rwl != sign_j");
        rwl = 4'b1010;   // slot0=0, slot1=1, slot2=0, slot3=1
        #1;
        for (int s=0; s<K_MAX_T; s++) begin
            logic exp_s;
            logic [IC_BITS_T-1:0] exp_j;
            exp_s = expected_xnor_bit(exp_sign[s], rwl[s]);
            exp_j = expected_xnor_J(exp_J[s], rwl[s]);
            check_equal_bit($sformatf("slot%d sign_eq", s), rbl_signs[s], exp_s);
            check_equal_vec($sformatf("slot%d J bits", s), get_rbl_J(s), exp_j);
        end

        // -----------------------------------------------------------------
        // 4. Precharge (forces all rbl to 1)
        // -----------------------------------------------------------------
        $display("Test 4: Precharge");
        precharge = 1;
        #1;
        for (int s=0; s<K_MAX_T; s++) begin
            check_equal_bit($sformatf("precharge sign_eq slot%d", s), rbl_signs[s], 1'b1);
            check_equal_vec($sformatf("precharge J bits slot%d", s), get_rbl_J(s), {IC_BITS_T{1'b1}});
        end
        precharge = 0;

        // -----------------------------------------------------------------
        // 5. Overwrite slot 1 and verify only that slot changes
        // -----------------------------------------------------------------
        $display("Test 5: Overwrite slot 1");
        write_slot_and_verify(1, 5, 4'b0011, 1'b1);
        // Verify all slots again with rwl=1 to ensure only slot1 changed
        rwl = 4'b1111;
        #1;
        for (int s=0; s<K_MAX_T; s++) begin
            logic exp_s;
            logic [IC_BITS_T-1:0] exp_j;
            exp_s = expected_xnor_bit(exp_sign[s], 1); // rwl=1 => stored
            exp_j = expected_xnor_J(exp_J[s], 1);
            check_equal_bit($sformatf("slot%d sign after overwrite", s), rbl_signs[s], exp_s);
            check_equal_vec($sformatf("slot%d J after overwrite", s), get_rbl_J(s), exp_j);
        end

        // -----------------------------------------------------------------
        // Summary
        // -----------------------------------------------------------------
        $display("=== Test Summary ===");
        $display("Total checks passed: %0d", test_pass_count);
        if (error_count == 0) begin
            $display("=== All tests PASSED ===");
        end else begin
            $error("=== %0d tests FAILED ===", error_count);
            $finish(1);
        end
        $finish;
    end

    // ─── Optional: hierarchical probe of internal q of sign bitcell for slot0 ───
    // This is for debugging only (not synthesizable)
    initial begin
        #100;
        // Probe q of sign bitcell slot0 (if simulation supports it)
        // You can add this to waveform viewer manually.
        // Path: tb_sparse_compute_memory_advanced.u_dut.gen_slot[0].u_sign.u_J.q
        $display("To check internal q of sign bitcell slot0, add:");
        $display("  tb_sparse_compute_memory_advanced.u_dut.gen_slot[0].u_sign.u_J.q");
    end

    // VCD dump
    initial begin
        $dumpfile("tb_sparse_compute_memory_advanced.vcd");
        $dumpvars(0, tb_sparse_compute_memory_advanced);
    end

endmodule
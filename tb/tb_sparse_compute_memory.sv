// tb_sparse_compute_memory.sv
// Self-checking testbench for sparse_compute_memory
`timescale 1ns/1ps

module tb_sparse_compute_memory;

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

    // Helper: extract a slot's J bits
    function automatic logic [IC_BITS_T-1:0] get_rbl_J(input int slot);
        return rbl_J[slot*IC_BITS_T +: IC_BITS_T];
    endfunction

    // Helper: compute expected XNOR of a stored vector with a single-bit rwl
    function automatic logic [IC_BITS_T-1:0] expected_xnor_J(
        input logic [IC_BITS_T-1:0] stored,
        input logic                 rwl
    );
        return (rwl ? stored : ~stored);
    endfunction

    // Helper: compute expected sign_eq for a slot (XNOR of stored sign_j and rwl)
    function automatic logic expected_sign_eq(input logic sign_j, input logic rwl);
        return ~(sign_j ^ rwl);
    endfunction

    // Check tasks
    task check_equal_vec(input string msg, input logic [IC_BITS_T-1:0] a, input logic [IC_BITS_T-1:0] b);
        if (a !== b) begin
            $error("%s: expected %h, got %h", msg, b, a);
            error_count++;
        end else begin
            $display("  %s: OK", msg);
        end
    endtask

    task check_equal_bit(input string msg, input logic a, input logic b);
        if (a !== b) begin
            $error("%s: expected %b, got %b", msg, b, a);
            error_count++;
        end else begin
            $display("  %s: OK", msg);
        end
    endtask

    // ─── Stored expected data (for reference) ──────────────────────────
    logic [COL_W_T-1:0]          exp_col [K_MAX_T];
    logic [IC_BITS_T-1:0]        exp_J  [K_MAX_T];
    logic                        exp_sign[K_MAX_T];

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
        $display("=== Starting sparse_compute_memory test ===");

        // -----------------------------------------------------------------
        // 1. Write all slots with known values
        // -----------------------------------------------------------------
        $display("Test 1: Write all slots");
        // Slot 0: col=2, J=5 (0101), sign=1
        @(posedge clk);
        wr_en = 1; wr_slot = 0; wr_col_idx = 2; wr_J = 5; wr_sign_j = 1;
        exp_col[0] = 2; exp_J[0] = 5; exp_sign[0] = 1;
        @(posedge clk);
        // Slot 1: col=7, J=10 (1010), sign=0
        wr_en = 1; wr_slot = 1; wr_col_idx = 7; wr_J = 10; wr_sign_j = 0;
        exp_col[1] = 7; exp_J[1] = 10; exp_sign[1] = 0;
        @(posedge clk);
        // Slot 2: col=3, J=7 (0111), sign=1
        wr_en = 1; wr_slot = 2; wr_col_idx = 3; wr_J = 7; wr_sign_j = 1;
        exp_col[2] = 3; exp_J[2] = 7; exp_sign[2] = 1;
        @(posedge clk);
        // Slot 3: col=0, J=15 (1111), sign=0
        wr_en = 1; wr_slot = 3; wr_col_idx = 0; wr_J = 15; wr_sign_j = 0;
        exp_col[3] = 0; exp_J[3] = 15; exp_sign[3] = 0;
        @(posedge clk);
        wr_en = 0;

        // -----------------------------------------------------------------
        // 2. Read col_idx (with correct timing)
        // -----------------------------------------------------------------
        $display("Test 2: Read col_idx");
        // Read slot 0
        @(posedge clk);
        rd_col_en = 1; rd_slot = 0;
        @(posedge clk);
        // Now rd_col_valid is high, rd_col_idx is valid
        check_equal_vec("Read col_idx slot0", rd_col_idx, exp_col[0]);
        if (rd_col_valid !== 1) $error("rd_col_valid not high for slot0");
        rd_col_en = 0;

        // Read slot 2
        @(posedge clk);
        rd_col_en = 1; rd_slot = 2;
        @(posedge clk);
        check_equal_vec("Read col_idx slot2", rd_col_idx, exp_col[2]);
        if (rd_col_valid !== 1) $error("rd_col_valid not high for slot2");
        rd_col_en = 0;

        // -----------------------------------------------------------------
        // 3. Compute: rwl == sign_j (all sign_eq should be 1)
        //    rwl = {slot3=0, slot2=1, slot1=0, slot0=1} → 4'b0101
        // -----------------------------------------------------------------
        $display("Test 3: Compute with rwl == sign_j");
        rwl = 4'b0101;   // slot0=1, slot1=0, slot2=1, slot3=0
        #1;
        for (int s=0; s<K_MAX_T; s++) begin
            // sign_eq
            logic exp_s = expected_sign_eq(exp_sign[s], rwl[s]);
            check_equal_bit($sformatf("slot%d sign_eq", s), rbl_signs[s], exp_s);
            // J bits
            logic [IC_BITS_T-1:0] exp_j = expected_xnor_J(exp_J[s], rwl[s]);
            check_equal_vec($sformatf("slot%d J bits", s), get_rbl_J(s), exp_j);
        end

        // -----------------------------------------------------------------
        // 4. Compute: rwl != sign_j (all sign_eq should be 0)
        //    rwl = {slot3=1, slot2=0, slot1=1, slot0=0} → 4'b1010
        // -----------------------------------------------------------------
        $display("Test 4: Compute with rwl != sign_j");
        rwl = 4'b1010;   // slot0=0, slot1=1, slot2=0, slot3=1
        #1;
        for (int s=0; s<K_MAX_T; s++) begin
            logic exp_s = expected_sign_eq(exp_sign[s], rwl[s]);
            check_equal_bit($sformatf("slot%d sign_eq", s), rbl_signs[s], exp_s);
            logic [IC_BITS_T-1:0] exp_j = expected_xnor_J(exp_J[s], rwl[s]);
            check_equal_vec($sformatf("slot%d J bits", s), get_rbl_J(s), exp_j);
        end

        // -----------------------------------------------------------------
        // 5. Precharge (forces all rbl to 1)
        // -----------------------------------------------------------------
        $display("Test 5: Precharge");
        precharge = 1;
        #1;
        for (int s=0; s<K_MAX_T; s++) begin
            check_equal_bit($sformatf("precharge sign_eq slot%d", s), rbl_signs[s], 1'b1);
            check_equal_vec($sformatf("precharge J bits slot%d", s), get_rbl_J(s), {IC_BITS_T{1'b1}});
        end
        precharge = 0;

        // -----------------------------------------------------------------
        // 6. Write slot 1 new value, then re‑read and compute
        // -----------------------------------------------------------------
        $display("Test 6: Overwrite slot 1");
        @(posedge clk);
        wr_en = 1; wr_slot = 1; wr_col_idx = 5; wr_J = 3; wr_sign_j = 1;
        exp_col[1] = 5; exp_J[1] = 3; exp_sign[1] = 1;
        @(posedge clk);
        wr_en = 0;

        // Read col_idx slot1
        @(posedge clk);
        rd_col_en = 1; rd_slot = 1;
        @(posedge clk);
        check_equal_vec("Read col_idx slot1 after update", rd_col_idx, exp_col[1]);
        if (rd_col_valid !== 1) $error("rd_col_valid not high for slot1 after update");
        rd_col_en = 0;

        // Compute with rwl == new sign_j (slot1=1, others don't matter)
        rwl = 4'b0010;   // only slot1=1
        #1;
        for (int s=0; s<K_MAX_T; s++) begin
            logic exp_s = expected_sign_eq(exp_sign[s], rwl[s]);
            check_equal_bit($sformatf("slot%d sign_eq after update", s), rbl_signs[s], exp_s);
            logic [IC_BITS_T-1:0] exp_j = expected_xnor_J(exp_J[s], rwl[s]);
            check_equal_vec($sformatf("slot%d J bits after update", s), get_rbl_J(s), exp_j);
        end

        // -----------------------------------------------------------------
        // Summary
        // -----------------------------------------------------------------
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
        $dumpfile("tb_sparse_compute_memory.vcd");
        $dumpvars(0, tb_sparse_compute_memory);
    end

endmodule
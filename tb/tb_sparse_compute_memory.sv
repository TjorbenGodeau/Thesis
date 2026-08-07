// tb_sparse_compute_memory.sv
`timescale 1ns/1ps

module tb_sparse_compute_memory;

    import dsb_pkg::*;
    import dsb_sparse_pkg::*;

    // Parameters for test (K_MAX=4, IC_BITS=4, COL_W=3)
    localparam int K_MAX_T  = 4;
    localparam int IC_BITS_T = 4;
    localparam int COL_W_T   = 3;  // up to 8 oscillators

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

    // Helper to read a slot's J bits and sign from rbl arrays
    function automatic logic [IC_BITS_T-1:0] get_rbl_J(input int slot);
        return rbl_J[slot*IC_BITS_T +: IC_BITS_T];
    endfunction

    function automatic logic get_rbl_sign(input int slot);
        return rbl_signs[slot];
    endfunction

    // Test control
    int error_count = 0;

    task check_equal(input string msg, input logic a, input logic b);
        if (a !== b) begin
            $error("%s: expected %b, got %b", msg, b, a);
            error_count++;
        end else begin
            $display("  %s: OK", msg);
        end
    endtask

    task check_equal_vec(input string msg, input logic [IC_BITS_T-1:0] a, input logic [IC_BITS_T-1:0] b);
        if (a !== b) begin
            $error("%s: expected %b, got %b", msg, b, a);
            error_count++;
        end else begin
            $display("  %s: OK", msg);
        end
    endtask

    // Test sequence
    initial begin
        // Initialize
        wr_en      = 0;
        wr_slot    = 0;
        wr_col_idx = 0;
        wr_J       = 0;
        wr_sign_j  = 0;
        rd_col_en  = 0;
        rd_slot    = 0;
        precharge  = 0;
        rwl        = 0;

        repeat (2) @(posedge clk);
        $display("=== Starting sparse_compute_memory test ===");

        // --------------------------------------------------------------------
        // 1. Write all slots with known values
        // --------------------------------------------------------------------
        $display("Test 1: Write slots");
        // Slot 0: col=2, J=5 (0101), sign=1
        @(posedge clk);
        wr_en = 1; wr_slot = 0; wr_col_idx = 2; wr_J = 5; wr_sign_j = 1;
        @(posedge clk);
        // Slot 1: col=7, J=10 (1010), sign=0
        wr_en = 1; wr_slot = 1; wr_col_idx = 7; wr_J = 10; wr_sign_j = 0;
        @(posedge clk);
        // Slot 2: col=3, J=7 (0111), sign=1
        wr_en = 1; wr_slot = 2; wr_col_idx = 3; wr_J = 7; wr_sign_j = 1;
        @(posedge clk);
        // Slot 3: col=0, J=15 (1111), sign=0
        wr_en = 1; wr_slot = 3; wr_col_idx = 0; wr_J = 15; wr_sign_j = 0;
        @(posedge clk);
        wr_en = 0;

        // --------------------------------------------------------------------
        // 2. Read col_idx (registered output, 1-cycle latency)
        // --------------------------------------------------------------------
        $display("Test 2: Read col_idx");
        // Read slot 0
        @(posedge clk);
        rd_col_en = 1; rd_slot = 0;
        @(posedge clk);
        check_equal("Read col_idx slot0", rd_col_idx, 2);
        if (rd_col_valid !== 1) $error("rd_col_valid not high");
        rd_col_en = 0;

        // Read slot 2
        @(posedge clk);
        rd_col_en = 1; rd_slot = 2;
        @(posedge clk);
        check_equal("Read col_idx slot2", rd_col_idx, 3);
        if (rd_col_valid !== 1) $error("rd_col_valid not high");
        rd_col_en = 0;

        // Read slot not written (slot 0 after another write? Not needed)

        // --------------------------------------------------------------------
        // 3. Compute with rwl matching stored sign_j
        //    rwl = [1,0,1,0]  (slot0:1, slot1:0, slot2:1, slot3:0)
        //    All rbl_signs should be 1, rbl_J should equal stored J.
        // --------------------------------------------------------------------
        $display("Test 3: Compute with rwl == sign_j (all sign_eq=1)");
        @(posedge clk);
        precharge = 0;
        rwl = 4'b1010;  // slot3=0, slot2=1, slot1=0, slot0=1
        #1; // allow combinational settling
        // Check each slot
        check_equal("slot0 sign_eq", get_rbl_sign(0), 1);
        check_equal_vec("slot0 J bits", get_rbl_J(0), 5);
        check_equal("slot1 sign_eq", get_rbl_sign(1), 1);
        check_equal_vec("slot1 J bits", get_rbl_J(1), 10);
        check_equal("slot2 sign_eq", get_rbl_sign(2), 1);
        check_equal_vec("slot2 J bits", get_rbl_J(2), 7);
        check_equal("slot3 sign_eq", get_rbl_sign(3), 1);
        check_equal_vec("slot3 J bits", get_rbl_J(3), 15);

        // --------------------------------------------------------------------
        // 4. Compute with rwl != sign_j for some slots
        //    rwl = [0,1,0,1]  (slot0:0 mismatch, slot1:1 mismatch, slot2:0 mismatch, slot3:1 mismatch)
        //    sign_eq should be 0 for all slots.
        //    J bits: XNOR(J_bit, rwl)
        //      slot0: J=0101, rwl=0 -> ~0101 = 1010
        //      slot1: J=1010, rwl=1 -> 1010 (identity)
        //      slot2: J=0111, rwl=0 -> ~0111 = 1000
        //      slot3: J=1111, rwl=1 -> 1111
        // --------------------------------------------------------------------
        $display("Test 4: Compute with rwl != sign_j (all sign_eq=0)");
        @(posedge clk);
        rwl = 4'b0101;  // slot3=1, slot2=0, slot1=1, slot0=0
        #1;
        check_equal("slot0 sign_eq", get_rbl_sign(0), 0);
        check_equal_vec("slot0 J bits (XNOR with 0)", get_rbl_J(0), ~5); // 1010
        check_equal("slot1 sign_eq", get_rbl_sign(1), 0);
        check_equal_vec("slot1 J bits (XNOR with 1)", get_rbl_J(1), 10); // 1010
        check_equal("slot2 sign_eq", get_rbl_sign(2), 0);
        check_equal_vec("slot2 J bits (XNOR with 0)", get_rbl_J(2), ~7); // 1000
        check_equal("slot3 sign_eq", get_rbl_sign(3), 0);
        check_equal_vec("slot3 J bits (XNOR with 1)", get_rbl_J(3), 15); // 1111

        // --------------------------------------------------------------------
        // 5. Precharge forces all rbl to 1
        // --------------------------------------------------------------------
        $display("Test 5: Precharge");
        @(posedge clk);
        precharge = 1;
        #1;
        for (int k = 0; k < K_MAX_T; k++) begin
            check_equal($sformatf("precharge sign_eq slot%d", k), get_rbl_sign(k), 1);
            check_equal_vec($sformatf("precharge J bits slot%d", k), get_rbl_J(k), {IC_BITS_T{1'b1}});
        end
        precharge = 0;

        // --------------------------------------------------------------------
        // 6. Write to a slot while others unchanged, then re-read to verify
        //    Also test that rwl doesn't affect write.
        // --------------------------------------------------------------------
        $display("Test 6: Write slot 1 new value");
        @(posedge clk);
        wr_en = 1; wr_slot = 1; wr_col_idx = 5; wr_J = 3; wr_sign_j = 1;
        @(posedge clk);
        wr_en = 0;
        // Read col_idx slot1
        @(posedge clk);
        rd_col_en = 1; rd_slot = 1;
        @(posedge clk);
        check_equal("Read col_idx slot1 after update", rd_col_idx, 5);
        rd_col_en = 0;
        @(posedge clk);
        // Compute with rwl == new sign_j (rwl=1 for slot1)
        rwl = 4'b0000; // set all 0 for simplicity, but set slot1=1
        rwl[1] = 1;
        #1;
        check_equal("slot1 sign_eq after update", get_rbl_sign(1), 1);
        check_equal_vec("slot1 J bits after update", get_rbl_J(1), 3);

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

    // VCD dump for waveform inspection
    initial begin
        $dumpfile("tb_sparse_compute_memory.vcd");
        $dumpvars(0, tb_sparse_compute_memory);
    end

endmodule
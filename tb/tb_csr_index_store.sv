// tb_csr_index_store.sv
`timescale 1ns/1ps

module tb_csr_index_store;

    import dsb_pkg::*;
    import dsb_sparse_pkg::*;

    // -------------------------------------------------------------------------
    // Local parameters for testing (small, to keep sim fast)
    // -------------------------------------------------------------------------
    localparam int N_TEST          = 4;
    localparam int MAX_NNZ_TEST    = 16;
    localparam int COL_IDX_W_TEST  = $clog2(N_TEST);   // 2
    localparam int MAIN_ADDR_W_TEST = 4;
    localparam int CSR_ADDR_W_TEST = $clog2(MAX_NNZ_TEST); // 4
    localparam int CSR_ENTRY_W_TEST = COL_IDX_W_TEST + MAIN_ADDR_W_TEST; // 6
    localparam int ROW_COUNT_W_TEST = $clog2(MAX_NNZ_TEST); // 4
    localparam int ROW_PTR_W_TEST   = CSR_ADDR_W_TEST + ROW_COUNT_W_TEST; // 8

    // Clock
    logic clk = 0;
    always #5 clk = ~clk;

    // DUT signals
    logic                              rst_n;
    logic                              rp_wr_en;
    logic [COL_IDX_W_TEST-1:0]         rp_wr_row;
    logic [ROW_PTR_W_TEST-1:0]         rp_wr_data;
    logic                              rp_rd_en;
    logic [COL_IDX_W_TEST-1:0]         rp_rd_row;
    logic [ROW_PTR_W_TEST-1:0]         rp_rd_data;
    logic                              rp_rd_valid;
    logic                              en_wr_en;
    logic [CSR_ADDR_W_TEST-1:0]        en_wr_addr;
    logic [CSR_ENTRY_W_TEST-1:0]       en_wr_data;
    logic                              en_rd_en;
    logic [CSR_ADDR_W_TEST-1:0]        en_rd_addr;
    logic [CSR_ENTRY_W_TEST-1:0]       en_rd_data;
    logic                              en_rd_valid;

    // DUT
    csr_index_store #(
        .N_P          (N_TEST),
        .MAX_NNZ_P    (MAX_NNZ_TEST),
        .COL_IDX_W_P  (COL_IDX_W_TEST),
        .MAIN_ADDR_WP (MAIN_ADDR_W_TEST),
        .CSR_ADDR_W_P (CSR_ADDR_W_TEST),
        .CSR_ENTRY_WP (CSR_ENTRY_W_TEST),
        .ROW_PTR_W_P  (ROW_PTR_W_TEST),
        .ROW_COUNT_WP (ROW_COUNT_W_TEST)
    ) u_dut (
        .clk         (clk),
        .rst_n       (rst_n),
        .rp_wr_en    (rp_wr_en),
        .rp_wr_row   (rp_wr_row),
        .rp_wr_data  (rp_wr_data),
        .rp_rd_en    (rp_rd_en),
        .rp_rd_row   (rp_rd_row),
        .rp_rd_data  (rp_rd_data),
        .rp_rd_valid (rp_rd_valid),
        .en_wr_en    (en_wr_en),
        .en_wr_addr  (en_wr_addr),
        .en_wr_data  (en_wr_data),
        .en_rd_en    (en_rd_en),
        .en_rd_addr  (en_rd_addr),
        .en_rd_data  (en_rd_data),
        .en_rd_valid (en_rd_valid)
    );

    // Test control
    int error_count = 0;

    task check_equal_val(input string msg, input logic a, input logic b);
        if (a !== b) begin
            $error("%s: expected %b, got %b", msg, b, a);
            error_count++;
        end else begin
            $display("  %s: OK", msg);
        end
    endtask

    task check_equal_vec(input string msg, input logic [31:0] a, input logic [31:0] b);
        if (a !== b) begin
            $error("%s: expected 0x%h, got 0x%h", msg, b, a);
            error_count++;
        end else begin
            $display("  %s: OK", msg);
        end
    endtask

    // ─── Row pointer write ────────────────────────────────────────────────
    task automatic write_rp(input int row, input logic [ROW_PTR_W_TEST-1:0] data);
        @(posedge clk);
        rp_wr_en = 1; rp_wr_row = row; rp_wr_data = data;
        @(posedge clk);
        rp_wr_en = 0;
    endtask

    // ─── Row pointer read (correct timing) ──────────────────────────────
    task automatic read_rp(input int row, input logic [ROW_PTR_W_TEST-1:0] expected);
        @(posedge clk);
        rp_rd_en = 1; rp_rd_row = row;
        @(posedge clk);
        // Now rd_data and rd_valid are valid (latched on this edge)
        check_equal_vec($sformatf("RP read row %0d", row), rp_rd_data, expected);
        check_equal_val($sformatf("RP valid row %0d", row), rp_rd_valid, 1);
        rp_rd_en = 0;
    endtask

    // ─── Entry write ─────────────────────────────────────────────────────
    task automatic write_en(input int addr, input logic [CSR_ENTRY_W_TEST-1:0] data);
        @(posedge clk);
        en_wr_en = 1; en_wr_addr = addr; en_wr_data = data;
        @(posedge clk);
        en_wr_en = 0;
    endtask

    // ─── Entry read (correct timing) ────────────────────────────────────
    task automatic read_en(input int addr, input logic [CSR_ENTRY_W_TEST-1:0] expected);
        @(posedge clk);
        en_rd_en = 1; en_rd_addr = addr;
        @(posedge clk);
        check_equal_vec($sformatf("EN read addr %0d", addr), en_rd_data, expected);
        check_equal_val($sformatf("EN valid addr %0d", addr), en_rd_valid, 1);
        en_rd_en = 0;
    endtask

    // ─── Test sequence ──────────────────────────────────────────────────
    initial begin
        // Initialize
        rst_n      = 0;
        rp_wr_en   = 0;
        rp_rd_en   = 0;
        en_wr_en   = 0;
        en_rd_en   = 0;
        repeat (2) @(posedge clk);
        rst_n      = 1;
        repeat (2) @(posedge clk);
        $display("=== Starting csr_index_store test ===");

        // -----------------------------------------------------------------
        // 1. Write all row pointers
        // -----------------------------------------------------------------
        $display("Test 1: Write all row pointers");
        for (int i = 0; i < N_TEST; i++) begin
            automatic logic [ROW_PTR_W_TEST-1:0] rp_data;
            rp_data = {CSR_ADDR_W_TEST'(i*2), ROW_COUNT_W_TEST'(i+1)};
            write_rp(i, rp_data);
            $display("  Wrote row %0d: start=%0d, count=%0d", i, i*2, i+1);
        end

        // -----------------------------------------------------------------
        // 2. Read them back and verify
        // -----------------------------------------------------------------
        $display("Test 2: Read row pointers");
        for (int i = 0; i < N_TEST; i++) begin
            automatic logic [ROW_PTR_W_TEST-1:0] expected;
            expected = {CSR_ADDR_W_TEST'(i*2), ROW_COUNT_W_TEST'(i+1)};
            read_rp(i, expected);
        end

        // -----------------------------------------------------------------
        // 3. Write entries
        // -----------------------------------------------------------------
        $display("Test 3: Write entries");
        for (int a = 0; a < 8; a++) begin
            automatic logic [CSR_ENTRY_W_TEST-1:0] en_data;
            en_data = {COL_IDX_W_TEST'(a+1), MAIN_ADDR_W_TEST'(a*3)};
            write_en(a, en_data);
            $display("  Wrote EN addr %0d: col=%0d, maddr=%0d", a, a+1, a*3);
        end

        // -----------------------------------------------------------------
        // 4. Read entries back
        // -----------------------------------------------------------------
        $display("Test 4: Read entries");
        for (int a = 0; a < 8; a++) begin
            automatic logic [CSR_ENTRY_W_TEST-1:0] expected;
            expected = {COL_IDX_W_TEST'(a+1), MAIN_ADDR_W_TEST'(a*3)};
            read_en(a, expected);
        end

        // -----------------------------------------------------------------
        // 5. Overwrite entry 3 and 5
        // -----------------------------------------------------------------
        $display("Test 5: Overwrite entry 3 and 5");
        write_en(3, {COL_IDX_W_TEST'(7), MAIN_ADDR_W_TEST'(42)});
        write_en(5, {COL_IDX_W_TEST'(2), MAIN_ADDR_W_TEST'(99)});
        read_en(3, {COL_IDX_W_TEST'(7), MAIN_ADDR_W_TEST'(42)});
        read_en(5, {COL_IDX_W_TEST'(2), MAIN_ADDR_W_TEST'(99)});
        read_en(2, {COL_IDX_W_TEST'(3), MAIN_ADDR_W_TEST'(6)});
        read_en(4, {COL_IDX_W_TEST'(5), MAIN_ADDR_W_TEST'(12)});

        // -----------------------------------------------------------------
        // 6. Overwrite row pointer 2
        // -----------------------------------------------------------------
        $display("Test 6: Overwrite row pointer 2");
        write_rp(2, {CSR_ADDR_W_TEST'(9), ROW_COUNT_W_TEST'(5)});
        read_rp(2, {CSR_ADDR_W_TEST'(9), ROW_COUNT_W_TEST'(5)});
        read_rp(1, {CSR_ADDR_W_TEST'(2), ROW_COUNT_W_TEST'(2)});
        read_rp(3, {CSR_ADDR_W_TEST'(6), ROW_COUNT_W_TEST'(4)});

        // -----------------------------------------------------------------
        // 7. Reset behaviour
        // -----------------------------------------------------------------
        $display("Test 7: Reset behaviour");
        read_rp(0, {CSR_ADDR_W_TEST'(0), ROW_COUNT_W_TEST'(1)});
        @(posedge clk);
        rst_n = 0;
        repeat (2) @(posedge clk);
        check_equal_val("RP valid after reset", rp_rd_valid, 0);
        check_equal_vec("RP data after reset", rp_rd_data, 0);
        check_equal_val("EN valid after reset", en_rd_valid, 0);
        check_equal_vec("EN data after reset", en_rd_data, 0);
        @(posedge clk);
        rst_n = 1;
        repeat (2) @(posedge clk);
        read_rp(0, {CSR_ADDR_W_TEST'(0), ROW_COUNT_W_TEST'(1)});
        read_en(0, {COL_IDX_W_TEST'(1), MAIN_ADDR_W_TEST'(0)});

        // -----------------------------------------------------------------
        // 8. Read uninitialised entry (validity only)
        // -----------------------------------------------------------------
        $display("Test 8: Read uninitialised entry (validity only)");
        @(posedge clk);
        en_rd_en = 1; en_rd_addr = 15;
        @(posedge clk);
        check_equal_val("EN valid for uninitialised", en_rd_valid, 1);
        en_rd_en = 0;

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
        $dumpfile("tb_csr_index_store.vcd");
        $dumpvars(0, tb_csr_index_store);
    end

endmodule
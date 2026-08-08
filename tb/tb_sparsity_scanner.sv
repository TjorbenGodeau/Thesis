// tb_sparsity_scanner.sv
`timescale 1ns/1ps

module tb_sparsity_scanner;

    import dsb_pkg::*;
    import dsb_sparse_pkg::*;

    // -------------------------------------------------------------------------
    // Local parameters for a tiny test system
    // -------------------------------------------------------------------------
    localparam int N_TEST          = 4;
    localparam int MAX_EDGES_TEST  = 8;
    localparam int MAX_NNZ_TEST    = 2 * MAX_EDGES_TEST;
    localparam int COL_IDX_W_TEST  = $clog2(N_TEST);             // 2
    localparam int MAIN_ADDR_W_TEST = $clog2(MAX_EDGES_TEST + 2*N_TEST); // 4
    localparam int CSR_ADDR_W_TEST = $clog2(MAX_NNZ_TEST + 1);   // 4
    localparam int CSR_ENTRY_W_TEST = COL_IDX_W_TEST + MAIN_ADDR_W_TEST; // 6
    localparam int ROW_COUNT_W_TEST = $clog2(MAX_NNZ_TEST + 1);  // 4
    localparam int ROW_PTR_W_TEST   = CSR_ADDR_W_TEST + ROW_COUNT_W_TEST; // 8
    localparam int MAIN_WORD_W_TEST = ( (COL_IDX_W_TEST+COL_IDX_W_TEST+IC_BITS) > XY_W ) ?
                                       (COL_IDX_W_TEST+COL_IDX_W_TEST+IC_BITS) : XY_W;
    localparam int MAIN_MEM_WORDS_TEST = MAX_EDGES_TEST + 2*N_TEST;

    // Clock
    logic clk = 0;
    always #5 clk = ~clk;

    // Reset
    logic rst_n = 1;

    // Scanner signals
    logic                     scan_start;
    logic [CSR_ADDR_W_TEST-1:0] num_edges;
    logic [COL_IDX_W_TEST-1:0]  active_n;
    logic                     scan_done;
    logic [CSR_ADDR_W_TEST-1:0] nnz_total;

    // Main memory signals
    logic                     mm_wr_en;
    logic [MAIN_ADDR_W_TEST-1:0] mm_wr_addr;
    logic [MAIN_WORD_W_TEST-1:0] mm_wr_data;
    logic                     mm_rd_en;
    logic [MAIN_ADDR_W_TEST-1:0] mm_rd_addr;
    logic [MAIN_WORD_W_TEST-1:0] mm_rd_data;
    logic                     mm_rd_valid;

    // CSR store write ports
    logic                     rp_wr_en;
    logic [COL_IDX_W_TEST-1:0] rp_wr_row;
    logic [ROW_PTR_W_TEST-1:0] rp_wr_data;
    logic                     en_wr_en;
    logic [CSR_ADDR_W_TEST-1:0] en_wr_addr;
    logic [CSR_ENTRY_W_TEST-1:0] en_wr_data;

    // CSR store read ports (testbench reads)
    logic                     rp_rd_en;
    logic [COL_IDX_W_TEST-1:0] rp_rd_row;
    logic [ROW_PTR_W_TEST-1:0] rp_rd_data;
    logic                     rp_rd_valid;
    logic                     en_rd_en;
    logic [CSR_ADDR_W_TEST-1:0] en_rd_addr;
    logic [CSR_ENTRY_W_TEST-1:0] en_rd_data;
    logic                     en_rd_valid;

    // -------------------------------------------------------------------------
    // Instantiations
    // -------------------------------------------------------------------------
    main_memory #(
        .WORDS_P  (MAIN_MEM_WORDS_TEST),
        .WORD_W_P (MAIN_WORD_W_TEST),
        .ADDR_W_P (MAIN_ADDR_W_TEST)
    ) u_mm (
        .clk     (clk),
        .rst_n   (rst_n),
        .wr_en   (mm_wr_en),
        .wr_addr (mm_wr_addr),
        .wr_data (mm_wr_data),
        .rd_en   (mm_rd_en),
        .rd_addr (mm_rd_addr),
        .rd_data (mm_rd_data),
        .rd_valid(mm_rd_valid)
    );

    csr_index_store #(
        .N_P          (N_TEST),
        .MAX_NNZ_P    (MAX_NNZ_TEST),
        .COL_IDX_W_P  (COL_IDX_W_TEST),
        .MAIN_ADDR_WP (MAIN_ADDR_W_TEST),
        .CSR_ADDR_W_P (CSR_ADDR_W_TEST),
        .CSR_ENTRY_WP (CSR_ENTRY_W_TEST),
        .ROW_PTR_W_P  (ROW_PTR_W_TEST),
        .ROW_COUNT_WP (ROW_COUNT_W_TEST)
    ) u_csr (
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

    sparsity_scanner #(
        .N_P          (N_TEST),
        .IC_BITS_P    (IC_BITS),
        .MAIN_ADDR_WP (MAIN_ADDR_W_TEST),
        .MAIN_WORD_WP (MAIN_WORD_W_TEST),
        .COL_IDX_W_P  (COL_IDX_W_TEST),
        .CSR_ADDR_W_P (CSR_ADDR_W_TEST),
        .CSR_ENTRY_WP (CSR_ENTRY_W_TEST),
        .ROW_PTR_W_P  (ROW_PTR_W_TEST),
        .ROW_COUNT_WP (ROW_COUNT_W_TEST)
    ) u_scan (
        .clk        (clk),
        .rst_n      (rst_n),
        .start      (scan_start),
        .num_edges  (num_edges),
        .active_n   (active_n),
        .mm_rd_en   (mm_rd_en),
        .mm_rd_addr (mm_rd_addr),
        .mm_rd_data (mm_rd_data),
        .mm_rd_valid(mm_rd_valid),
        .rp_wr_en   (rp_wr_en),
        .rp_wr_row  (rp_wr_row),
        .rp_wr_data (rp_wr_data),
        .en_wr_en   (en_wr_en),
        .en_wr_addr (en_wr_addr),
        .en_wr_data (en_wr_data),
        .done       (scan_done),
        .nnz_total  (nnz_total)
    );

    // -------------------------------------------------------------------------
    // Testbench helpers
    // -------------------------------------------------------------------------
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

    // Pack edge word: {row, col, weight}
    // word = {row[COL_IDX_W-1:0], col[COL_IDX_W-1:0], weight[IC_BITS-1:0]}
    // Pad with zeros at MSB if MAIN_WORD_W is wider.
    function automatic logic [MAIN_WORD_W_TEST-1:0] pack_edge(
        input int r, c, w
    );
        logic [MAIN_WORD_W_TEST-1:0] p;
        p = '0;
        // Place row at top bits
        p[MAIN_WORD_W_TEST-1 -: COL_IDX_W_TEST] = COL_IDX_W_TEST'(r);
        // Place col just below row
        p[IC_BITS_P + COL_IDX_W_TEST -: COL_IDX_W_TEST] = COL_IDX_W_TEST'(c);
        // Place weight at LSB
        p[IC_BITS_P-1:0] = IC_BITS'(signed'(w));
        return p;
    endfunction

    // Write to main memory
    task mm_write(input int addr, input logic [MAIN_WORD_W_TEST-1:0] data);
        @(posedge clk);
        mm_wr_en = 1; mm_wr_addr = addr; mm_wr_data = data;
        @(posedge clk);
        mm_wr_en = 0;
    endtask

    // ─── Read CSR row pointer (correct timing) ──────────────────────────────
    task automatic check_rp(input int row, input int expected_start, input int expected_count);
        logic [ROW_PTR_W_TEST-1:0] expected;
        expected = {CSR_ADDR_W_TEST'(expected_start), ROW_COUNT_W_TEST'(expected_count)};
        @(posedge clk);
        rp_rd_en = 1; rp_rd_row = row;
        @(posedge clk);
        // Now rd_data and rd_valid are valid
        check_equal_vec($sformatf("rp_rd_data row %0d", row), rp_rd_data, expected);
        check_equal_val($sformatf("rp_rd_valid row %0d", row), rp_rd_valid, 1);
        rp_rd_en = 0;
    endtask

    // ─── Read CSR entry (correct timing) ────────────────────────────────────
    task automatic check_en(input int addr, input int expected_col, input int expected_maddr);
        logic [CSR_ENTRY_W_TEST-1:0] expected;
        expected = {COL_IDX_W_TEST'(expected_col), MAIN_ADDR_W_TEST'(expected_maddr)};
        @(posedge clk);
        en_rd_en = 1; en_rd_addr = addr;
        @(posedge clk);
        check_equal_vec($sformatf("en_rd_data addr %0d", addr), en_rd_data, expected);
        check_equal_val($sformatf("en_rd_valid addr %0d", addr), en_rd_valid, 1);
        en_rd_en = 0;
    endtask

    // -------------------------------------------------------------------------
    // Test cases
    // -------------------------------------------------------------------------
    initial begin
        // ---- All declarations must be at the top ----
        int EDGES[0:2][2];
        int num_edges_test;
        int active_n_test;
        int e;   // loop variable

        // Initialize signals
        mm_wr_en   = 0;
        mm_wr_addr = 0;
        mm_wr_data = 0;
        scan_start = 0;
        num_edges  = 0;
        active_n   = 0;
        rp_rd_en   = 0;
        en_rd_en   = 0;

        // Initialize test data
        EDGES = '{{0,1}, {0,2}, {1,3}};
        num_edges_test = 3;
        active_n_test = 4;

        repeat (4) @(posedge clk);
        rst_n = 1;
        repeat (2) @(posedge clk);
        $display("=== Starting sparsity_scanner test ===");

        // ── Test 1: Normal graph (3 edges, N=4) ──────────────────────────────
        $display("Test 1: Normal graph with 3 edges");

        // Preload main memory
        for (e = 0; e < num_edges_test; e++) begin
            mm_write(MAIN_EDGE_BASE + e,
                     pack_edge(EDGES[e][0], EDGES[e][1], 1));
        end

        // Set scanner inputs
        num_edges = CSR_ADDR_W_TEST'(num_edges_test);
        active_n  = COL_IDX_W_TEST'(active_n_test);

        // Start scanner
        @(posedge clk);
        scan_start = 1;
        @(posedge clk);
        scan_start = 0;

        // Wait for done (with timeout)
        repeat (1000) @(posedge clk) if (scan_done) break;
        if (!scan_done) $fatal("Timeout waiting for scan_done");

        $display("scan_done, nnz_total = %0d", nnz_total);

        // Verify row pointers
        check_rp(0, 0, 2);   // row0: entries 0,1
        check_rp(1, 2, 2);   // row1: entries 2,3
        check_rp(2, 4, 1);   // row2: entry 4
        check_rp(3, 5, 1);   // row3: entry 5

        // Verify entries
        // Expected order: edges processed in given order.
        // Edge0 (0,1): forward 0->1 at addr 0, reverse 1->0 at addr 2
        // Edge1 (0,2): forward 0->2 at addr 1, reverse 2->0 at addr 4
        // Edge2 (1,3): forward 1->3 at addr 3, reverse 3->1 at addr 5
        check_en(0, 1, MAIN_EDGE_BASE+0);
        check_en(1, 2, MAIN_EDGE_BASE+1);
        check_en(2, 0, MAIN_EDGE_BASE+0);
        check_en(3, 3, MAIN_EDGE_BASE+2);
        check_en(4, 0, MAIN_EDGE_BASE+1);
        check_en(5, 1, MAIN_EDGE_BASE+2);

        // --------------------------------------------------------------------
        // Test 2: Empty edge list
        // --------------------------------------------------------------------
        $display("Test 2: Empty edge list (num_edges=0)");
        // Set num_edges=0, active_n=4
        num_edges = 0;
        active_n  = COL_IDX_W_TEST'(4);
        @(posedge clk);
        scan_start = 1;
        @(posedge clk);
        scan_start = 0;
        repeat (1000) @(posedge clk) if (scan_done) break;
        if (!scan_done) $fatal("Timeout waiting for empty scan");

        // All rows should have count=0, start=0
        for (int i = 0; i < active_n_test; i++) begin
            check_rp(i, 0, 0);
        end
        // nnz_total should be 0
        check_equal_val("nnz_total", nnz_total, 0);

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
        $dumpfile("tb_sparsity_scanner.vcd");
        $dumpvars(0, tb_sparsity_scanner);
    end

endmodule
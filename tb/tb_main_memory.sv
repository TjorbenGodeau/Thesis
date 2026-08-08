// tb_main_memory.sv
`timescale 1ns/1ps

module tb_main_memory;

    import dsb_pkg::*;
    import dsb_sparse_pkg::*;

    // Parameters (use default from dsb_sparse_pkg)
    localparam WORDS_P  = MAIN_MEM_WORDS;
    localparam WORD_W_P = MAIN_WORD_W;
    localparam ADDR_W_P = MAIN_ADDR_W;

    // Clock
    logic clk = 0;
    always #5 clk = ~clk;

    // DUT signals
    logic                  wr_en;
    logic [ADDR_W_P-1:0]   wr_addr;
    logic [WORD_W_P-1:0]   wr_data;

    logic                  rd_en;
    logic [ADDR_W_P-1:0]   rd_addr;
    logic [WORD_W_P-1:0]   rd_data;
    logic                  rd_valid;

    // DUT
    main_memory u_dut (
        .clk     (clk),
        .rst_n   (1'b1),   // no reset needed for memory content, but we tie high
        .wr_en   (wr_en),
        .wr_addr (wr_addr),
        .wr_data (wr_data),
        .rd_en   (rd_en),
        .rd_addr (rd_addr),
        .rd_data (rd_data),
        .rd_valid(rd_valid)
    );

    // Test control
    int error_count = 0;

    task check_equal(input string msg, input logic [WORD_W_P-1:0] a, input logic [WORD_W_P-1:0] b);
        if (a !== b) begin
            $error("%s: expected 0x%h, got 0x%h", msg, b, a);
            error_count++;
        end else begin
            $display("  %s: OK", msg);
        end
    endtask

    task check_equal_val(input string msg, input logic a, input logic b);
        if (a !== b) begin
            $error("%s: expected %b, got %b", msg, b, a);
            error_count++;
        end else begin
            $display("  %s: OK", msg);
        end
    endtask

    // Helper to write and read back (one address)
    task write_read(input logic [ADDR_W_P-1:0] addr,
                    input logic [WORD_W_P-1:0] data);
        // Write
        @(posedge clk);
        wr_en = 1; wr_addr = addr; wr_data = data;
        @(posedge clk);
        wr_en = 0;
        // Read back (issue read one cycle after write, wait for valid)
        @(posedge clk);
        rd_en = 1; rd_addr = addr;
        @(posedge clk);
        rd_en = 0;
        // Data appears on next cycle with valid
        @(posedge clk);
        check_equal($sformatf("Read back addr 0x%0h", addr), rd_data, data);
        check_equal_val($sformatf("rd_valid for addr 0x%0h", addr), rd_valid, 1);
    endtask

    initial begin
        // Initialize
        wr_en = 0;
        wr_addr = 0;
        wr_data = 0;
        rd_en = 0;
        rd_addr = 0;

        repeat (2) @(posedge clk);
        $display("=== Starting main_memory test ===");

        // ── Test 1: Write/read single address ──────────────────────────────
        $display("Test 1: Write/read single address");
        write_read(5, 32'hDEADBEEF);

        // ── Test 2: Write/read multiple addresses ──────────────────────────
        $display("Test 2: Multiple writes and reads");
        write_read(10, 32'hA5A5A5A5);
        write_read(20, 32'h5A5A5A5A);
        write_read(30, 32'hFFFFFFFF);

        // Read them again (no writes in between)
        @(posedge clk);
        rd_en = 1; rd_addr = 10;
        @(posedge clk);
        rd_en = 0;
        @(posedge clk);
        check_equal("Read addr 10 again", rd_data, 32'hA5A5A5A5);
        check_equal_val("rd_valid addr 10", rd_valid, 1);

        @(posedge clk);
        rd_en = 1; rd_addr = 20;
        @(posedge clk);
        rd_en = 0;
        @(posedge clk);
        check_equal("Read addr 20 again", rd_data, 32'h5A5A5A5A);
        check_equal_val("rd_valid addr 20", rd_valid, 1);

        @(posedge clk);
        rd_en = 1; rd_addr = 30;
        @(posedge clk);
        rd_en = 0;
        @(posedge clk);
        check_equal("Read addr 30 again", rd_data, 32'hFFFFFFFF);
        check_equal_val("rd_valid addr 30", rd_valid, 1);

        // ── Test 3: rd_valid only pulses when rd_en asserted ──────────────
        $display("Test 3: rd_valid only on read enable");
        // Issue a read
        @(posedge clk);
        rd_en = 1; rd_addr = 10;
        @(posedge clk);
        rd_en = 0;
        @(posedge clk);
        // rd_valid should be 1 here
        if (rd_valid !== 1) $error("rd_valid should be 1 after read enable");
        else $display("  rd_valid pulse OK");
        // Next cycle, without rd_en, rd_valid should be 0
        @(posedge clk);
        if (rd_valid !== 0) $error("rd_valid should be 0 without rd_en");
        else $display("  rd_valid deasserts OK");

        // ── Test 4: Overwrite and verify ──────────────────────────────────
        $display("Test 4: Overwrite address");
        write_read(5, 32'h12345678); // new data
        // Read again to confirm overwritten
        @(posedge clk);
        rd_en = 1; rd_addr = 5;
        @(posedge clk);
        rd_en = 0;
        @(posedge clk);
        check_equal("Overwritten addr 5", rd_data, 32'h12345678);
        check_equal_val("rd_valid overwritten", rd_valid, 1);

        // ── Test 5: Uninitialized read (should be X or 0, but we can't check) ──
        // But we can read an address never written. In simulation, uninitialized memory
        // is X, but some tools initialize to 0. We'll just check that rd_valid still pulses.
        $display("Test 5: Read uninitialized address (validity only)");
        @(posedge clk);
        rd_en = 1; rd_addr = 100; // never written
        @(posedge clk);
        rd_en = 0;
        @(posedge clk);
        if (rd_valid !== 1) $error("rd_valid should be 1 for uninitialized read");
        else $display("  rd_valid OK for uninitialized");
        // rd_data may be X, but we can't check.

        // ── Test 6: All zeros and all ones ──────────────────────────────────
        $display("Test 6: Write zeros and ones");
        write_read(200, 0);
        write_read(201, {WORD_W_P{1'b1}});

        // ── Test 7: Random values (just a few) ─────────────────────────────
        $display("Test 7: Random values");
        write_read(50, 32'hC0FFEE12);
        write_read(51, 32'hDEADBEEF);
        write_read(52, 32'h12345678);

        // ── Summary ──────────────────────────────────────────────────────────
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
        $dumpfile("tb_main_memory.vcd");
        $dumpvars(0, tb_main_memory);
    end

endmodule
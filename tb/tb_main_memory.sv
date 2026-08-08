// tb_main_memory.sv
// Self-checking testbench for main_memory with reduced size for simulation.
`timescale 1ns/1ps

module tb_main_memory;

    import dsb_pkg::*;
    import dsb_sparse_pkg::*;

    // ─── Override parameters for simulation (small size) ────────────────────
    localparam int TEST_WORDS = 64;          // small enough for waveform viewing
    localparam int TEST_WORD_W = 32;         // keep default word width
    localparam int TEST_ADDR_W = $clog2(TEST_WORDS);

    // Clock
    logic clk = 0;
    always #5 clk = ~clk;

    // DUT signals (use test parameters)
    logic                  wr_en;
    logic [TEST_ADDR_W-1:0] wr_addr;
    logic [TEST_WORD_W-1:0] wr_data;

    logic                  rd_en;
    logic [TEST_ADDR_W-1:0] rd_addr;
    logic [TEST_WORD_W-1:0] rd_data;
    logic                  rd_valid;

    // DUT with overridden parameters
    main_memory #(
        .WORDS_P  (TEST_WORDS),
        .WORD_W_P (TEST_WORD_W),
        .ADDR_W_P (TEST_ADDR_W)
    ) u_dut (
        .clk     (clk),
        .rst_n   (1'b1),
        .wr_en   (wr_en),
        .wr_addr (wr_addr),
        .wr_data (wr_data),
        .rd_en   (rd_en),
        .rd_addr (rd_addr),
        .rd_data (rd_data),
        .rd_valid(rd_valid)
    );

    // ─── Test control ────────────────────────────────────────────────────────
    int error_count = 0;

    task check_equal(input string msg, input logic [TEST_WORD_W-1:0] a, input logic [TEST_WORD_W-1:0] b);
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

    // ─── Write/read task (correct timing) ──────────────────────────────────
    // 1) Assert wr_en and data on a clock edge.
    // 2) Deassert wr_en on the next edge.
    // 3) For read: assert rd_en on a clock edge.
    // 4) On the following edge, rd_data and rd_valid become valid;
    //    check them immediately, then deassert rd_en.
    task write_read(input logic [TEST_ADDR_W-1:0] addr,
                    input logic [TEST_WORD_W-1:0] data);
        // Write phase
        @(posedge clk);
        wr_en = 1; wr_addr = addr; wr_data = data;
        @(posedge clk);
        wr_en = 0;

        // Read phase – correct timing: assert rd_en
        @(posedge clk);
        rd_en = 1; rd_addr = addr;
        @(posedge clk);
        // Now rd_data and rd_valid are valid (they changed on this edge)
        check_equal($sformatf("Read back addr 0x%0h", addr), rd_data, data);
        check_equal_val($sformatf("rd_valid for addr 0x%0h", addr), rd_valid, 1);
        // Deassert rd_en (can be done now; will take effect next cycle)
        rd_en = 0;
    endtask

    // ─── Test sequence ──────────────────────────────────────────────────────
    initial begin
        // Initialize
        wr_en = 0; wr_addr = 0; wr_data = 0;
        rd_en = 0; rd_addr = 0;

        repeat (2) @(posedge clk);
        $display("=== Starting main_memory test (reduced size) ===");

        // Test 1: Write/read single address
        $display("Test 1: Write/read single address");
        write_read(5, 32'hDEADBEEF);

        // Test 2: Multiple writes and reads
        $display("Test 2: Multiple writes and reads");
        write_read(10, 32'hA5A5A5A5);
        write_read(20, 32'h5A5A5A5A);
        write_read(30, 32'hFFFFFFFF);

        // Read them again (no writes in between)
        @(posedge clk);
        rd_en = 1; rd_addr = 10;
        @(posedge clk);
        check_equal("Read addr 10 again", rd_data, 32'hA5A5A5A5);
        check_equal_val("rd_valid addr 10", rd_valid, 1);
        rd_en = 0;

        @(posedge clk);
        rd_en = 1; rd_addr = 20;
        @(posedge clk);
        check_equal("Read addr 20 again", rd_data, 32'h5A5A5A5A);
        check_equal_val("rd_valid addr 20", rd_valid, 1);
        rd_en = 0;

        @(posedge clk);
        rd_en = 1; rd_addr = 30;
        @(posedge clk);
        check_equal("Read addr 30 again", rd_data, 32'hFFFFFFFF);
        check_equal_val("rd_valid addr 30", rd_valid, 1);
        rd_en = 0;

        // Test 3: rd_valid only pulses when rd_en is asserted
        $display("Test 3: rd_valid only on read enable");
        @(posedge clk);
        rd_en = 1; rd_addr = 10;
        @(posedge clk);
        // rd_valid should be 1 here
        if (rd_valid !== 1) $error("rd_valid should be 1 after read enable");
        else $display("  rd_valid pulse OK");
        rd_en = 0;
        // Next cycle, rd_valid should be 0
        @(posedge clk);
        if (rd_valid !== 0) $error("rd_valid should be 0 without rd_en");
        else $display("  rd_valid deasserts OK");

        // Test 4: Overwrite and verify
        $display("Test 4: Overwrite address");
        write_read(5, 32'h12345678);
        // Read again to confirm overwritten
        @(posedge clk);
        rd_en = 1; rd_addr = 5;
        @(posedge clk);
        check_equal("Overwritten addr 5", rd_data, 32'h12345678);
        check_equal_val("rd_valid overwritten", rd_valid, 1);
        rd_en = 0;

        // Test 5: Read uninitialized address (validity only)
        $display("Test 5: Read uninitialized address (validity only)");
        @(posedge clk);
        rd_en = 1; rd_addr = 100; // within 0..63? Actually 100 is out of range for 64 words.
        // Use an address within range but never written, e.g., 40.
        rd_addr = 40;
        @(posedge clk);
        if (rd_valid !== 1) $error("rd_valid should be 1 for uninitialized read");
        else $display("  rd_valid OK for uninitialized");
        rd_en = 0;
        // rd_data may be X, but we can't check.

        // Test 6: All zeros and all ones
        $display("Test 6: Write zeros and ones");
        write_read(50, 0);
        write_read(51, {TEST_WORD_W{1'b1}});

        // Test 7: Random values
        $display("Test 7: Random values");
        write_read(52, 32'hC0FFEE12);
        write_read(53, 32'hDEADBEEF);
        write_read(54, 32'h12345678);

        // Summary
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
        $dumpfile("tb_main_memory.vcd");
        $dumpvars(0, tb_main_memory);
    end

endmodule
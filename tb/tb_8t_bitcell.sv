// tb_8t_bitcell.sv
`timescale 1ns/1ps

module tb_8t_bitcell;

    // Parameters
    localparam CLK_PERIOD = 10; // ns

    // Signals
    logic clk;
    logic wwl;
    logic wbl;
    logic rwl;
    logic precharge;
    logic rbl_out;

    // Instantiate the bitcell
    dsb_8t_bitcell u_bitcell (
        .clk       (clk),
        .wwl       (wwl),
        .wbl       (wbl),
        .rwl       (rwl),
        .precharge (precharge),
        .rbl_out   (rbl_out)
    );

    // Clock generation
    initial begin
        clk = 0;
        forever #(CLK_PERIOD/2) clk = ~clk;
    end

    // Test sequence
    initial begin
        // Initialize signals
        wwl       = 0;
        wbl       = 0;
        rwl       = 0;
        precharge = 0;

        // Wait a few cycles
        repeat (2) @(posedge clk);
        $display("=== Starting 8T bitcell test ===");

        // ------------------------------------------------------------
        // Test 1: Write 0 and verify
        // ------------------------------------------------------------
        $display("Test 1: Write 0");
        @(posedge clk);
        wwl = 1;
        wbl = 0;
        @(posedge clk);
        wwl = 0;   // finish write

        // Read with rwl=0 (should get XNOR(0,0)=1)
        #1; // allow combinational settle
        rwl = 0;
        precharge = 0;
        #1;
        if (rbl_out !== 1) $error("Write 0, rwl=0: expected 1, got %b", rbl_out);
        else $display("  Write 0, rwl=0 -> rbl_out=1 (OK)");

        // Read with rwl=1 (should get XNOR(0,1)=0)
        rwl = 1;
        #1;
        if (rbl_out !== 0) $error("Write 0, rwl=1: expected 0, got %b", rbl_out);
        else $display("  Write 0, rwl=1 -> rbl_out=0 (OK)");

        // ------------------------------------------------------------
        // Test 2: Write 1 and verify
        // ------------------------------------------------------------
        $display("Test 2: Write 1");
        @(posedge clk);
        wwl = 1;
        wbl = 1;
        @(posedge clk);
        wwl = 0;

        // Read with rwl=0 (should get XNOR(1,0)=0)
        rwl = 0;
        #1;
        if (rbl_out !== 0) $error("Write 1, rwl=0: expected 0, got %b", rbl_out);
        else $display("  Write 1, rwl=0 -> rbl_out=0 (OK)");

        // Read with rwl=1 (should get XNOR(1,1)=1)
        rwl = 1;
        #1;
        if (rbl_out !== 1) $error("Write 1, rwl=1: expected 1, got %b", rbl_out);
        else $display("  Write 1, rwl=1 -> rbl_out=1 (OK)");

        // ------------------------------------------------------------
        // Test 3: Precharge override (should force rbl_out=1)
        // ------------------------------------------------------------
        $display("Test 3: Precharge override");
        precharge = 1;
        rwl = 0;   // regardless of stored value
        #1;
        if (rbl_out !== 1) $error("Precharge: expected 1, got %b", rbl_out);
        else $display("  Precharge -> rbl_out=1 (OK)");
        precharge = 0;

        // ------------------------------------------------------------
        // Test 4: wwl normal mode override (should force rbl_out=1)
        // ------------------------------------------------------------
        $display("Test 4: wwl normal mode override");
        wwl = 1;
        rwl = 0;
        #1;
        if (rbl_out !== 1) $error("wwl normal mode: expected 1, got %b", rbl_out);
        else $display("  wwl=1 -> rbl_out=1 (OK)");
        wwl = 0;

        // ------------------------------------------------------------
        // Summary
        // ------------------------------------------------------------
        $display("=== Test finished ===");
        $finish;
    end

    // VCD dump (optional, for waveform inspection)
    initial begin
        $dumpfile("tb_8t_bitcell.vcd");
        $dumpvars(0, tb_8t_bitcell);
    end

endmodule
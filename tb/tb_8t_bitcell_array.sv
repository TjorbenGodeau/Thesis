// tb_8t_bitcell_array.sv
`timescale 1ns/1ps

module tb_8t_bitcell_array;

    // --------------------------------------------------------------
    // Parameters for the array size
    // --------------------------------------------------------------
    localparam int SLOTS = 4;        // number of word lines (slots)
    localparam int BITS  = 4;        // number of bit lines per slot (J width)

    // Clock
    logic clk = 0;
    always #5 clk = ~clk;

    // --------------------------------------------------------------
    // Signals for the array of bitcells
    // --------------------------------------------------------------
    logic                  precharge;
    logic [SLOTS-1:0]      wwl;         // write word line per slot
    logic [BITS-1:0]       wbl;         // write bit lines (shared across slots)
    logic [SLOTS-1:0]      rwl;         // read word line per slot
    logic [SLOTS*BITS-1:0] rbl_out;     // read bit lines (one per bit per slot)

    // --------------------------------------------------------------
    // Instantiate the bitcell array
    // --------------------------------------------------------------
    for (genvar s = 0; s < SLOTS; s++) begin : gen_slot
        for (genvar b = 0; b < BITS; b++) begin : gen_bit
            dsb_8t_bitcell u_bit (
                .clk       (clk),
                .wwl       (wwl[s]),
                .wbl       (wbl[b]),
                .rwl       (rwl[s]),
                .precharge (precharge),
                .rbl_out   (rbl_out[s*BITS + b])
            );
        end
    end

    // --------------------------------------------------------------
    // Helper: extract a slot's rbl vector
    // --------------------------------------------------------------
    function automatic logic [BITS-1:0] get_rbl(input int slot);
        return rbl_out[slot*BITS +: BITS];
    endfunction

    // --------------------------------------------------------------
    // Test control
    // --------------------------------------------------------------
    int error_count = 0;

    task check_equal(input string msg, input logic [BITS-1:0] a, input logic [BITS-1:0] b);
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

    // Helper to write a full slot (all bits)
    task write_slot(input int slot, input logic [BITS-1:0] data);
        @(posedge clk);
        wwl[slot] = 1;
        wbl = data;
        @(posedge clk);
        wwl[slot] = 0;
    endtask

    // Helper to compute expected XNOR for a slot given stored data and rwl
    function automatic logic [BITS-1:0] expected_xnor(
        input logic [BITS-1:0] stored,
        input logic            rwl
    );
        // XNOR of each bit: (stored_bit == rwl)
        return {BITS{(rwl)}} ^ stored;   // if rwl=0, XNOR = ~stored; if rwl=1, XNOR = stored
        // Actually: XNOR(a,b) = ~(a ^ b). For bits: a ^ b = a ^ rwl (since b=rwl)
        // So XNOR = ~(a ^ rwl) = ~a if rwl=0, else a if rwl=1.
        // In Verilog bitwise: return (rwl ? stored : ~stored);
    endfunction

    // --------------------------------------------------------------
    // Test sequence
    // --------------------------------------------------------------
    initial begin
        // Initialize all signals
        precharge = 0;
        wwl = 0;
        wbl = 0;
        rwl = 0;

        repeat (2) @(posedge clk);
        $display("=== Advanced 8T Bitcell Array Test ===");
        $display("Array size: %0d slots × %0d bits", SLOTS, BITS);

        // ------------------------------------------------------------
        // 1. Write different values to each slot
        // ------------------------------------------------------------
        $display("Test 1: Write slots with patterns");
        write_slot(0, 4'b0101);   // 5
        write_slot(1, 4'b1010);   // 10
        write_slot(2, 4'b0111);   // 7
        write_slot(3, 4'b1111);   // 15

        // ------------------------------------------------------------
        // 2. Read with rwl = 1 (should return stored value exactly)
        // ------------------------------------------------------------
        $display("Test 2: Read with rwl=1 (identity)");
        rwl = 4'b1111;   // all slots enabled with rwl=1
        #1;  // allow combinational settle
        check_equal("slot0 rwl=1", get_rbl(0), 4'b0101);
        check_equal("slot1 rwl=1", get_rbl(1), 4'b1010);
        check_equal("slot2 rwl=1", get_rbl(2), 4'b0111);
        check_equal("slot3 rwl=1", get_rbl(3), 4'b1111);

        // ------------------------------------------------------------
        // 3. Read with rwl = 0 (should return bitwise NOT of stored)
        // ------------------------------------------------------------
        $display("Test 3: Read with rwl=0 (inversion)");
        rwl = 4'b0000;
        #1;
        check_equal("slot0 rwl=0", get_rbl(0), ~4'b0101);
        check_equal("slot1 rwl=0", get_rbl(1), ~4'b1010);
        check_equal("slot2 rwl=0", get_rbl(2), ~4'b0111);
        check_equal("slot3 rwl=0", get_rbl(3), ~4'b1111);

        // ------------------------------------------------------------
        // 4. Mixed rwl per slot (e.g., slot0=1, slot1=0, slot2=1, slot3=0)
        // ------------------------------------------------------------
        $display("Test 4: Mixed rwl = 4'b0101 (slot0=1, slot1=0, slot2=1, slot3=0)");
        rwl = 4'b0101;   // bit order: rwl[3] = slot3, [2]=slot2, [1]=slot1, [0]=slot0
        #1;
        check_equal("slot0 rwl=1", get_rbl(0), expected_xnor(4'b0101, 1));
        check_equal("slot1 rwl=0", get_rbl(1), expected_xnor(4'b1010, 0));
        check_equal("slot2 rwl=1", get_rbl(2), expected_xnor(4'b0111, 1));
        check_equal("slot3 rwl=0", get_rbl(3), expected_xnor(4'b1111, 0));

        // ------------------------------------------------------------
        // 5. Precharge forces all rbl = 1
        // ------------------------------------------------------------
        $display("Test 5: Precharge override");
        precharge = 1;
        #1;
        for (int s = 0; s < SLOTS; s++) begin
            check_equal($sformatf("precharge slot%d", s), get_rbl(s), {BITS{1'b1}});
        end
        precharge = 0;

        // ------------------------------------------------------------
        // 6. Write one slot while others retain old data
        // ------------------------------------------------------------
        $display("Test 6: Overwrite slot 1 with new value 3");
        write_slot(1, 4'b0011);
        // Read back all slots with rwl=1
        rwl = 4'b1111;
        #1;
        check_equal("slot0 after overwrite", get_rbl(0), 4'b0101);
        check_equal("slot1 new value",       get_rbl(1), 4'b0011);
        check_equal("slot2 unchanged",       get_rbl(2), 4'b0111);
        check_equal("slot3 unchanged",       get_rbl(3), 4'b1111);

        // ------------------------------------------------------------
        // 7. Test wwl normal mode override (wwl=1 forces rbl=1)
        //    Verify that even during a write, the read output is forced high.
        // ------------------------------------------------------------
        $display("Test 7: wwl normal mode override");
        // First, set rwl=0 so without wwl we'd get inverted data.
        rwl = 4'b0000;
        #1;
        // Without wwl, rbl should be ~stored.
        check_equal("slot0 before wwl override", get_rbl(0), ~4'b0101);
        // Now assert wwl for slot0 (should force rbl_out=1)
        wwl[0] = 1;
        #1;
        check_equal("slot0 with wwl=1", get_rbl(0), {BITS{1'b1}});
        wwl[0] = 0;
        #1;
        // After wwl goes low, rbl returns to normal.
        check_equal("slot0 after wwl deassert", get_rbl(0), ~4'b0101);

        // ------------------------------------------------------------
        // 8. Test simultaneous write and read (should not affect rbl)
        //    Write slot2 while reading slot0 with rwl=1.
        // ------------------------------------------------------------
        $display("Test 8: Simultaneous write/read");
        @(posedge clk);
        wwl[2] = 1; wbl = 4'b1000;   // write to slot2
        rwl = 4'b0001;               // read only slot0
        @(posedge clk);
        wwl[2] = 0;
        #1;
        // slot0 should still be 0101; slot2 should be updated but not read.
        check_equal("slot0 during write", get_rbl(0), 4'b0101);
        // Now read slot2 with rwl=1 to confirm write
        rwl = 4'b0100;   // slot2 only
        #1;
        check_equal("slot2 after simultaneous write", get_rbl(2), 4'b1000);

        // ------------------------------------------------------------
        // Summary
        // ------------------------------------------------------------
        if (error_count == 0) begin
            $display("=== All tests PASSED ===");
        end else begin
            $error("=== %0d tests FAILED ===", error_count);
            $finish(1);
        end
        $finish;
    end

    // --------------------------------------------------------------
    // VCD dump for waveform inspection
    // --------------------------------------------------------------
    initial begin
        $dumpfile("tb_8t_bitcell_array.vcd");
        $dumpvars(0, tb_8t_bitcell_array);
    end

endmodule
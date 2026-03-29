`timescale 1ns/1ps

module sram_1r1w_tb;

    // parameters for the DUT
    parameter int DEPTH = 16;
    parameter int WIDTH = 8;

    logic clk;
    logic [DEPTH-1:0] wwl;
    logic [WIDTH-1:0] wbl;
    logic [DEPTH-1:0] rwl;
    logic [WIDTH-1:0] rbl;

    // DUT instantiation
    sram_1r1w #(
        .DEPTH(DEPTH),
        .WIDTH(WIDTH)
    ) dut (
        .clk(clk),
        .wwl(wwl),
        .wbl(wbl),
        .rwl(rwl),
        .rbl(rbl)
    );

    // clock generator
    initial begin
        clk = 0;
        forever #5 clk = ~clk;
    end

    // one-hot helper
    function automatic logic [DEPTH-1:0] onehot(input int idx);
        logic [DEPTH-1:0] tmp;
        tmp = '0;
        if (idx >= 0 && idx < DEPTH)
            tmp[idx] = 1'b1;
        return tmp;
    endfunction

    task automatic write_word(input int row, input logic [WIDTH-1:0] data);
        begin
            wwl = onehot(row);
            wbl = data;
            @(posedge clk);
            wwl = '0;
            wbl = '0;
        end
    endtask

    task automatic read_word(input int row, output logic [WIDTH-1:0] data);
        begin
            rwl = onehot(row);
            @(posedge clk);
            data = rbl;
            rwl = '0;
        end
    endtask

    initial begin
        int i;
        logic [WIDTH-1:0] expected;
        logic [WIDTH-1:0] read_data;

        wwl = '0;
        wbl = '0;
        rwl = '0;

        // wait a couple of cycles for reset domain stability
        repeat (2) @(posedge clk);

        // fill SRAM with pattern 0xA5 + row index
        $display("[%0t] START WRITE PHASE", $time);
        for (i = 0; i < DEPTH; i++) begin
            expected = {WIDTH{1'b0}} | (8'hA5 + i);
            write_word(i, expected);
        end

        // empty read to capture release state
        @(posedge clk);
        if (rbl !== {WIDTH{1'b1}}) begin
            $error("Expected idle rbl = {{WIDTH{1'b1}}}, got %h", rbl);
        end

        // read back all locations
        $display("[%0t] START READ BACK PHASE", $time);
        for (i = 0; i < DEPTH; i++) begin
            expected = {WIDTH{1'b0}} | (8'hA5 + i);
            read_word(i, read_data);
            if (read_data !== expected) begin
                $error("Data mismatch at row %0d: expected %02h, got %02h", i, expected, read_data);
            end else begin
                $display("row %0d => %02h OK", i, read_data);
            end
        end

        // test concurrent write and read to different rows
        $display("[%0t] CONCURRENT R/W PHASE", $time);
        wwl = onehot(3);
        wbl = {WIDTH{1'b0}} + 8'hFF;
        rwl = onehot(5);
        @(posedge clk);
        wwl = '0;
        wbl = '0;
        rwl = '0;

        // next cycle read row 3 and row 5 values
        read_word(3, read_data);
        if (read_data !== 8'hFF) $error("Row 3 expected 0xFF, got %02h", read_data);

        read_word(5, read_data);
        expected = {WIDTH{1'b0}} | (8'hA5 + 5);
        if (read_data !== expected) $error("Row 5 expected %02h, got %02h", expected, read_data);

        $display("[%0t] TEST PASSED", $time);
        $finish;
    end

endmodule

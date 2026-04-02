module 8t_bitcell(
    input logic clk,        
    input logic wwl,        // write word line (normal mode)
    input logic wbl,        // write bit line
    input logic rwl,        // read word line (comput mode: = sign(x_i))
    input logic precharge,
    output logic rbl_out    // XNOR(stored_bit, rwl)
);

    logic q, qb;

    // --- Write port ---
    always_ff @(posedge clk) begin
        if (wwl) begin
            q <= wbl;
            qb <= ~wbl;
        end
    end

    // --- Read / Compute port (combinational discharge model) ---
    // RBL discharges when XNOR(q, rwl) = 1 (i.e. q == rwl)
    always_comb begin
        if (precharge) rbl_out = 1'b1;
        else if (wwl) rbl_out = 1'b1;   // normal mode: no compute discharge
        else rbl_out = ~(q ^ rwl);       // XNOR
    end

endmodule

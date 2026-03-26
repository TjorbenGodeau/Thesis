module sram_1r1w #(
    parameter int DEPTH  = 256,
    parameter int WIDTH  = 32
)(
    input  logic                clk,

    // write port
    input  logic [DEPTH-1:0]    wwl,
    input  logic [WIDTH-1:0]    wbl,

    // read port
    input  logic [DEPTH-1:0]    rwl,
    output logic [WIDTH-1:0]    rbl
);

logic [WIDTH-1:0] mem [DEPTH];
logic [WIDTH-1:0] read_q;

//
// helper: convert one-hot to index
//
function int onehot2index(logic [DEPTH-1:0] vec);
    for (int i=0; i<DEPTH; i++)
        if (vec[i]) return i;
    return 0;
endfunction

//
// WRITE
//
always_ff @(posedge clk) begin
    if (|wwl) begin
        int wrow = onehot2index(wwl);
        mem[wrow] <= wbl;
    end
end

//
// READ
//
always_ff @(posedge clk) begin
    if (|rwl) begin
        int rrow = onehot2index(rwl);
        read_q <= mem[rrow];
    end
end

//
// Precharged bitline abstraction
//
always_comb rbl = (|rwl) ? read_q : '1;

//
// Sanity checks
//
always_comb begin
    assert ($onehot0(wwl)) else $error("Multiple write wordlines active");
    assert ($onehot0(rwl)) else $error("Multiple read wordlines active");
end

endmodule
// filepath: c:/Users/tjorb/Documents/Thesis/rtl/UpdateMomentum.sv
// UpdateMomentum: Compute y_next = y - dt*(p*x - Jx)
// Implements Eq(3) from the GbSB model

module UpdateMomentum #(
    parameter int N = 8,
    parameter int DATA_WIDTH = 32,
    parameter int FRAC_WIDTH = 16
) (
    input  logic clk,
    input  logic rst_n,
    input  logic valid_in,
    
    input  logic [N-1:0][DATA_WIDTH-1:0] y,
    input  logic [N-1:0][DATA_WIDTH-1:0] x,
    input  logic [N-1:0][DATA_WIDTH-1:0] Jx,
    input  logic [N-1:0][DATA_WIDTH-1:0] p,
    input  logic [DATA_WIDTH-1:0] dt,
    
    output logic [N-1:0][DATA_WIDTH-1:0] y_next,
    output logic valid_out
);

    logic [N-1:0][DATA_WIDTH-1:0] px;
    logic [N-1:0][DATA_WIDTH-1:0] px_minus_Jx;
    logic [N-1:0][DATA_WIDTH-1:0] dt_term;
    logic [N-1:0][DATA_WIDTH-1:0] y_next_comb;
    logic valid_out_comb;

    // Stage 1: Compute p*x for each oscillator
    generate
        for (genvar i = 0; i < N; i++) begin : gen_px
            FixedPointMult #(
                .DATA_WIDTH(DATA_WIDTH),
                .FRAC_WIDTH(FRAC_WIDTH)
            ) mult_px (
                .a(p[i]),
                .b(x[i]),
                .product(px[i])
            );
        end
    endgenerate

    // Stage 2: Compute (p*x - Jx) for each oscillator
    generate
        for (genvar i = 0; i < N; i++) begin : gen_diff
            FixedPointSub #(
                .DATA_WIDTH(DATA_WIDTH),
                .FRAC_WIDTH(FRAC_WIDTH)
            ) sub_diff (
                .a(px[i]),
                .b(Jx[i]),
                .difference(px_minus_Jx[i])
            );
        end
    endgenerate

    // Stage 3: Compute dt * (p*x - Jx) for each oscillator
    generate
        for (genvar i = 0; i < N; i++) begin : gen_dt_term
            FixedPointMult #(
                .DATA_WIDTH(DATA_WIDTH),
                .FRAC_WIDTH(FRAC_WIDTH)
            ) mult_dt (
                .a(dt),
                .b(px_minus_Jx[i]),
                .product(dt_term[i])
            );
        end
    endgenerate

    // Stage 4: Compute y_next = y - dt*(p*x - Jx)
    generate
        for (genvar i = 0; i < N; i++) begin : gen_y_next
            FixedPointSub #(
                .DATA_WIDTH(DATA_WIDTH),
                .FRAC_WIDTH(FRAC_WIDTH)
            ) sub_y (
                .a(y[i]),
                .b(dt_term[i]),
                .difference(y_next_comb[i])
            );
        end
    endgenerate

    // Combinational output
    assign y_next = y_next_comb;
    assign valid_out_comb = valid_in;

    // Register for pipelining if needed
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            valid_out <= 1'b0;
        end else begin
            valid_out <= valid_out_comb;
        end
    end

endmodule
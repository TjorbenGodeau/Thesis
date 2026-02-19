// filepath: c:/Users/tjorb/Documents/Thesis/rtl/MatrixVectorMult.sv
// Matrix-Vector Multiplication: result = J * x
// Pipelined architecture for fixed-point computation

module MatrixVectorMult #(
    parameter int N = 8,           // Matrix dimension
    parameter int DATA_WIDTH = 32, // Fixed-point data width
    parameter int FRAC_WIDTH = 16  // Fractional bits
) (
    input  logic clk,
    input  logic rst_n,
    input  logic valid_in,
    
    // Inputs
    input  logic [N-1:0][N-1:0][DATA_WIDTH-1:0] J,
    input  logic [N-1:0][DATA_WIDTH-1:0] x,
    
    // Output
    output logic [N-1:0][DATA_WIDTH-1:0] result,
    output logic valid_out
);

    // Pipeline stages: accumulate products for each row
    logic [N-1:0][N-1:0][2*DATA_WIDTH-1:0] products;
    logic [N-1:0][2*DATA_WIDTH-1:0] accum;
    logic [N-1:0][DATA_WIDTH-1:0] result_next;
    logic valid_pipe;

    // Stage 1: Compute all J[i][j] * x[j] products
    always_comb begin
        for (int i = 0; i < N; i++) begin
            for (int j = 0; j < N; j++) begin
                // Fixed-point multiply: both operands have FRAC_WIDTH fractional bits
                products[i][j] = J[i][j] * x[j];
            end
        end
    end

    // Stage 2: Sum products for each row with fixed-point scaling
    always_comb begin
        for (int i = 0; i < N; i++) begin
            accum[i] = '0;
            for (int j = 0; j < N; j++) begin
                accum[i] += products[i][j];
            end
            // Scale back: remove extra FRAC_WIDTH bits from multiplication
            result_next[i] = accum[i][2*DATA_WIDTH-1-FRAC_WIDTH : DATA_WIDTH-FRAC_WIDTH];
        end
    end

    // Output register with valid flag
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            result <= '0;
            valid_out <= 1'b0;
        end else begin
            result <= result_next;
            valid_out <= valid_in;
        end
    end

endmodule
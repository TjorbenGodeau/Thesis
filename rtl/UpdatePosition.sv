// filepath: c:/Users/tjorb/Documents/Thesis/rtl/UpdatePosition.sv
// UpdatePosition: Update position x_next = x + dt*y
// Implements Eq (4) from GbSB model

module UpdatePosition #(
    parameter int N = 8,           // Number of oscillators
    parameter int DATA_WIDTH = 32, // Fixed-point data width
    parameter int FRAC_WIDTH = 16  // Fractional bits
) (
    input  logic clk,
    input  logic rst_n,
    input  logic valid_in,
    
    // Current state
    input  logic [N-1:0][DATA_WIDTH-1:0] x,
    input  logic [N-1:0][DATA_WIDTH-1:0] y,
    
    // Time step
    input  logic [DATA_WIDTH-1:0] dt,
    
    // Output
    output logic [N-1:0][DATA_WIDTH-1:0] x_next,
    output logic valid_out
);

    logic [N-1:0][DATA_WIDTH-1:0] x_next_comb;
    logic valid_out_ff;
    
    // Combinational computation of dt*y for each oscillator
    generate
        for (genvar i = 0; i < N; i++) begin : gen_update_x
            // Fixed-point multiplication: (dt * y[i]) >> FRAC_WIDTH
            // Then add to x[i]
            logic [2*DATA_WIDTH-1:0] product;
            
            always_comb begin
                product = dt * y[i];
                // Shift right by FRAC_WIDTH to extract integer result
                x_next_comb[i] = x[i] + product[2*DATA_WIDTH-1:FRAC_WIDTH];
            end
        end
    endgenerate
    
    // Pipeline register for valid signal
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            valid_out_ff <= 1'b0;
            x_next <= '0;
        end else begin
            valid_out_ff <= valid_in;
            x_next <= x_next_comb;
        end
    end
    
    assign valid_out = valid_out_ff;

endmodule
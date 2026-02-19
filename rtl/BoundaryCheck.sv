// filepath: c:/Users/tjorb/Documents/Thesis/rtl/BoundaryCheck.sv
// BoundaryCheck: Enforce inelastic wall collisions at |x| > 1
// Sets x = sign(x) and y = 0 for boundary violations

module BoundaryCheck #(
    parameter int N = 8,
    parameter int DATA_WIDTH = 32,
    parameter int FRAC_WIDTH = 16
) (
    input  logic clk,
    input  logic rst_n,
    input  logic valid_in,
    
    input  logic [N-1:0][DATA_WIDTH-1:0] x,
    input  logic [N-1:0][DATA_WIDTH-1:0] y,
    
    output logic [N-1:0][DATA_WIDTH-1:0] x_enforced,
    output logic [N-1:0][DATA_WIDTH-1:0] y_enforced,
    output logic [N-1:0] collision_mask,
    output logic any_collision,
    
    output logic valid_out
);
    
    logic [N-1:0][DATA_WIDTH-1:0] x_next, y_next;
    logic [N-1:0] collide;
    logic any_collide;
    
    // Fixed-point representation: 1.0 = 2^FRAC_WIDTH
    localparam logic [DATA_WIDTH-1:0] ONE = {{(DATA_WIDTH-FRAC_WIDTH-1){1'b0}}, 1'b1, {FRAC_WIDTH{1'b0}}};
    localparam logic [DATA_WIDTH-1:0] NEG_ONE = -ONE;
    
    always_comb begin
        any_collide = 1'b0;
        
        for (integer i = 0; i < N; i++) begin
            // Check if |x[i]| > 1.0
            logic [DATA_WIDTH-1:0] abs_x;
            abs_x = (x[i][DATA_WIDTH-1] == 1'b1) ? -x[i] : x[i];
            
            if (abs_x > ONE) begin
                collide[i] = 1'b1;
                any_collide = 1'b1;
                // Set x = sign(x)
                x_next[i] = (x[i][DATA_WIDTH-1] == 1'b1) ? NEG_ONE : ONE;
                // Set y = 0
                y_next[i] = '0;
            end else begin
                collide[i] = 1'b0;
                x_next[i] = x[i];
                y_next[i] = y[i];
            end
        end
    end
    
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            x_enforced <= '0;
            y_enforced <= '0;
            collision_mask <= '0;
            any_collision <= 1'b0;
            valid_out <= 1'b0;
        end else begin
            x_enforced <= x_next;
            y_enforced <= y_next;
            collision_mask <= collide;
            any_collision <= any_collide;
            valid_out <= valid_in;
        end
    end

endmodule
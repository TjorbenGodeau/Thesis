// dotprod_phase2_sparse.sv — includes sign_xi multiplication
`timescale 1ns/1ps
module dotprod_phase2_sparse
    import dsb_pkg::*;
    import dsb_sparse_pkg::*;
#(
    parameter int unsigned K_MAX_P    = K_MAX,
    parameter int unsigned IC_BITS_P  = IC_BITS,
    parameter int unsigned ACCUM_W_P  = ACCUM_W,
    parameter int unsigned K_CNT_W    = $clog2(K_MAX + 1)
)(
    input  logic                         clk,
    input  logic                         start,
    output logic                         done,
    input  logic                         accumulate,
    input  logic                         last_chunk,
    input  logic [K_CNT_W-1:0]           valid_count,
    input  logic                         sign_xi,
    input  logic [K_MAX_P*IC_BITS_P-1:0] xnor_J,
    input  logic [K_MAX_P-1:0]           sign_eq,
    output logic signed [ACCUM_W_P-1:0]  Jx_i
);

    logic signed [ACCUM_W_P-1:0] Jx_accum;

    always_ff @(posedge clk) begin
        done <= 1'b0;
        if (start) begin
            logic signed [ACCUM_W_P-1:0] chunk_sum;
            chunk_sum = '0;
            for (int k = 0; k < K_MAX_P; k++) begin
                if (k < int'(valid_count)) begin
                    logic [IC_BITS_P-1:0] xnor_word;
                    logic signed [IC_BITS_P-1:0] J_val;
                    logic signed [ACCUM_W_P-1:0] term;

                    xnor_word = xnor_J[k*IC_BITS_P +: IC_BITS_P];
                    // Recover J * sign_xi from xnor
                    J_val = sign_xi ? signed'(xnor_word) : signed'(~xnor_word);
                    // Compute J * σ_j = (sign_eq ? J : -J) * sign_xi
                    term = (sign_eq[k] ? J_val : -J_val);
                    if (~sign_xi) term = -term;   // multiply by -1 if sign_xi is -1
                    chunk_sum += ACCUM_W_P'(term);
                end
            end
            Jx_accum <= (accumulate ? Jx_accum : '0) + chunk_sum;
            done <= 1'b1;
        end
    end

    assign Jx_i = Jx_accum;

endmodule
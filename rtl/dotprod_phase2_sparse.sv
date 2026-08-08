// =============================================================================
// dotprod_phase2_sparse.sv – Correct per-term SACHI Equation (4)
// =============================================================================
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
    input  logic                         last_chunk,      // not used in correct model
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

                    // Reconstruct J from xnor and sign_xi (per SACHI Eq.4)
                    if (sign_xi)
                        J_val = signed'(xnor_word);      // xnor == J when sign_xi=+1
                    else
                        J_val = signed'(~xnor_word);     // xnor == ~J when sign_xi=-1

                    // The product J * sign_j * sign_xi = sign_eq ? J_val : -J_val
                    // because sign_eq = 1 iff sign_j == sign_xi.
                    term = sign_eq[k] ? J_val : -J_val;
                    chunk_sum += ACCUM_W_P'(term);
                end
            end

            // Accumulate: if accumulate=1, add to previous; else start fresh.
            Jx_accum <= (accumulate ? Jx_accum : '0) + chunk_sum;
            done <= 1'b1;
        end
    end

    assign Jx_i = Jx_accum;

endmodule
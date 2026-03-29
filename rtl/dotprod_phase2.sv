module dotprod_phase2
    import dsb_pkg::*;
    #(
        parameter int unsigned N_P = N,
        parameter int unsigned IC_BITS_P = IC_BITS,
        parameter int unsigned ACCUM_W_P = ACCUM_W
    )(
        input logic clk,
        input logic start,
        input logic sign_xi,
        input logic [N_P*IC_BITS_P-1:0] xnor_J,
        input logic [N_P-1:0] sign_eq,
        output logic signed [ACCUM_W_P-1:0] Jx_i,
        output logic done
    );

    always_ff @(posedge clk) begin
        done <= 1'b0;
        if (start) begin
            automatic logic signed [ACCUM_W_P-1:0] acc = '0;

            for (int j = 0; j < N_P; j++) begin
                automatic logic [IC_BITS_P-1:0] xnor_word = xnor_J[j*IC_BITS_P +: IC_BITS_P];
                // sign-equality gate
                automatic logic [IC_BITS_P-1:0] corrected = dign_eq[j] ? xnor_word : ~xnor_word;
                // sign-extend and axxumulate
                acc += ACCUM_W_P'(signed'(corrected));
            end

            // Global +1 per-term correction when sign(x_i) is negative
            if (!sign_xi) acc += ACCUM_W_P'(N_P);

            Jx_i <= acc;
            done <= 1'b1;
        end
    end
endmodule
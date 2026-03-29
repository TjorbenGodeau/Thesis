module schedule
    import dsb_pkg::*;
    #(
        parameter int unsigned A_BITS_P = A_BITS,
        parameter int unsigned STEP_W_P = STEP_W
    )(
        input logic [STEP_W_P-1:0] m,
        input logic [STEP_W_P-1:0] Nstep,
        input logic [A_BITS_P-1:0] a0_fp,
        output logic [A_BITS_P-1:0] a_m
    );

    // a(m) = a0_fp * (Nstep - m) / Nstep
    // Use extended precision to avoid overflow before the division
    logic [A_BITS_P+STEP_W_P-1:0] numerator;

    always_comb begin
        if (m >= Nstep) begin
            a_m = '0;
        end else begin
            numerator = a0_fp * (Nstep-m);
            a_m = A_BITS_P'(numerator / Nstep);
        end
    end

endmodule
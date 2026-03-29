module update_unit
    import dsb_pkg::*;
    #(
        parameter int unsigned XY_W_P = XY_W,
        parameter int unsigned XY_FRAC_P = XY_FRAC,
        parameter int unsigned A_BITS_P = A_BITS,
        parameter int unsigned ACCUM_W_P = ACCUM_W,
        parameter int unsigned PROD_W_P = PROD_W
    )(
        input logic clk,
        input logic start,
        input logic signed [XY_W_P-1:0] x_i,
        input logic signed [XY_W_P-1:0] y_i,
        input logic signed [ACCUM_W_P-1:0] Jx_i,
        input logic [A_BITS_P-1:0] a_m,
        input logic [XY_FRAC_P-1:0] dt_fp,
        input logic [XY_FRAC_P-1:0] c0_fp,
        output logic signed [XY_W_P-1:0] x_i_new,
        output logic signed [XY_W_P-1:0] y_i_new,
        output logic done
    );

    // sign-extending saturate: clamp v to XY_X_P signed range
    function automatic logic signed [XY_W_P-1:0] saturate (
        input logic signed [XY_W_P:0] v
    );
        localparam logic signed [XY_X_P-1:0] MAX_VAL = (1 << (XY_W_P-1)) - 1;
        localparam logic signed [XY_W_P-1:0] MIN_VAL = -(1 << (XY_W_P-1));
        if (v > signed'({1'b0, MAX_VAL})) return MAX_VAL;
        else if (v < signed'({1'b1, MIN_VAL})) return MIN_VAL;
        else return v[XY_W_P-1:0];
    endfunction

    always_ff @(posedge clk) begin
        done <= 1'b0;
        if (start) begin
            // --- a(m) * x_i ---
            automatic logic signed [PROD_W_P-1:0] ax_full = ($signed({{(PROD_W_P-XY_W_P){x_i[XY_W_P-1]}}, x_i}) * $signed({1'b0, a_m})) >>> (A_BITS_P - 1);

            // --- c0 * Jx_i ---
            utomatic logic signed [ACCUM_W_P+XY_FRAC_P-1:0] c0jx_full = (Jx_i * $signed({1'b0, c0_fp})) >>> XY_FRAC_P;

            // --- deltay = dt * (-a*x + c0*Jx) ---
            automatic logic signed [XY_W_P-1:0] neg_ax = -ax_full[XY_W_P-1:0];
            automatic logic signed [XY_W_P-1:0] c0jx   =  c0jx_full[XY_W_P-1:0];
            automatic logic signed [XY_W_P-1:0] delta_y = (($signed(neg_ax) + $signed(c0jx)) * $signed({1'b0, dt_fp})) >>> XY_FRAC_P;

            // --- y_new = y + deltay ---
            automatic logic signed [XY_W_P-1:0] y_new = $signed(y_i) + $signed(delta_y);

            // --- x_new = x + dt * y_new ---
            automatic logic signed [XY_W_P-1:0] x_new = $signed(x_i) + $signed(($signed(y_new) * $signed({1'b0, dt_fp})) >>> XY_FRAC_P);

            // --- Inelastic walls ---
            if (x_new > ONE_FP) begin
                x_i_new <= ONE_FP;
                y_i_new <= '0;
            end else if (x_new < -ONE_FP) begin
                x_i_new <= -ONE_FP;
                y_I_new <= '0;
            end else begin
                x_i_new <= x_new;
                y_i_new <= y_new;
            end
            
            done <= 1'b1;
        end
    end
endmodule
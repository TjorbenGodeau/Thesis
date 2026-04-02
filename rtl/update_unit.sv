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

    // ── Stage A: a(m) · x_i ───────────────────────────────────────────────────
    // Widen x_i to PROD_W_P bits by sign-extension before multiplying.
    // a_m is unsigned in [0,1], zero-extended by one bit to make it signed-safe.
    // Right-shift by (A_BITS_P-1) to renormalise back to Q(XY_INT).(XY_FRAC_P).
    logic signed [PROD_W_P-1:0]   ax_full;
    logic signed [XY_W_P-1:0]     ax_shifted;
    
    always_comb begin
        ax_full   = $signed({{(PROD_W_P-XY_W_P){x_i[XY_W_P-1]}}, x_i})
                    * $signed({1'b0, a_m});
        ax_shifted = XY_W_P'(ax_full >>> (A_BITS_P - 1));
    end
    
    // ── Stage B: c0 · Jx_i ────────────────────────────────────────────────────
    // Jx_i is a signed integer count of J-unit contributions (ACCUM_W_P wide).
    // c0_fp is Q0.XY_FRAC_P unsigned.  Product right-shifted by XY_FRAC_P
    // brings result back to Q(XY_INT).(XY_FRAC_P).
    logic signed [ACCUM_W_P+XY_FRAC_P-1:0] c0jx_full;
    logic signed [XY_W_P-1:0]               c0jx;
    
    always_comb begin
        c0jx_full = Jx_i * $signed({1'b0, c0_fp});
        c0jx      = XY_W_P'(c0jx_full >>> XY_FRAC_P);
    end
    
    // ── Stage C: Δy = dt · (−a·x + c0·Jx) ───────────────────────────────────
    // neg_ax negates the a·x term (Python: −a * x_i).
    // force  = −a·x + c0·Jx  is the net drive on momentum.
    // Multiplying by dt_fp and shifting gives Δy in Q(XY_INT).(XY_FRAC_P).
    logic signed [XY_W_P-1:0]     neg_ax;
    logic signed [XY_W_P-1:0]     force;
    logic signed [2*XY_W_P-1:0]   dy_full;
    logic signed [XY_W_P-1:0]     delta_y;
    
    always_comb begin
        neg_ax  = -ax_shifted;
        force   = neg_ax + c0jx;
        dy_full = $signed(force) * $signed({1'b0, dt_fp});
        delta_y = XY_W_P'(dy_full >>> XY_FRAC_P);
    end
    
    // ── Stage D: y_pre = y_i + Δy ─────────────────────────────────────────────
    // Direct addition — no saturation needed here because the inelastic wall
    // resets y to 0 whenever x hits ±1, keeping y bounded in practice.
    logic signed [XY_W_P-1:0] y_pre;
    
    always_comb begin
        y_pre = y_i + delta_y;
    end
    
    // ── Stage E: x_pre = x_i + dt · y_pre ────────────────────────────────────
    // Use the just-updated y_pre (Python: x += dt * y after y has been updated).
    logic signed [2*XY_W_P-1:0] dx_full;
    logic signed [XY_W_P-1:0]   x_pre;
    
    always_comb begin
        dx_full = $signed(y_pre) * $signed({1'b0, dt_fp});
        x_pre   = x_i + XY_W_P'(dx_full >>> XY_FRAC_P);
    end
    
    // ── Stage F: inelastic wall clamp ─────────────────────────────────────────
    // If |x_pre| > 1.0 (ONE_FP in fixed-point):
    //   x → sign(x_pre) × 1.0,  y → 0
    // ONE_FP is declared in dsb_pkg as  (1 << XY_FRAC).
    logic wall_hit_pos;
    logic wall_hit_neg;
    logic wall_hit;
    logic signed [XY_W_P-1:0] x_clamped;
    logic signed [XY_W_P-1:0] y_clamped;
    
    always_comb begin
        wall_hit_pos = (x_pre >  $signed(ONE_FP));
        wall_hit_neg = (x_pre < -$signed(ONE_FP));
        wall_hit     = wall_hit_pos | wall_hit_neg;
    
        unique case ({wall_hit_pos, wall_hit_neg})
        2'b10   : x_clamped =  ONE_FP;          // hit positive wall → +1.0
        2'b01   : x_clamped = -ONE_FP;          // hit negative wall → −1.0
        default : x_clamped =  x_pre;           // no wall hit → pass through
        endcase
    
        y_clamped = wall_hit ? '0 : y_pre;
    end
    
    // ── Output register — sample combinational result on start ─────────────────
    // All arithmetic above is live every cycle.  The always_ff simply captures
    // x_clamped and y_clamped into the output registers when start fires.
    always_ff @(posedge clk) begin
        done <= 1'b0;
        if (start) begin
        x_i_new <= x_clamped;
        y_i_new <= y_clamped;
        done    <= 1'b1;
        end
    end
endmodule
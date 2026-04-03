// =============================================================================
// 6.  Update Unit — momentum + position integration + inelastic wall
//
//     All intermediate computations are explicit named logic signals driven
//     by always_comb, forming real combinational hardware.  The always_ff
//     block only samples the final comb results into output registers when
//     start is asserted.
//
//     Datapath stages (all combinational, all active every cycle):
//
//       Stage A — a(m) · x_i
//         ax_full   [PROD_W_P]            = sign_ext(x_i) × zero_ext(a_m)
//         ax_shifted[XY_W_P]              = ax_full >>> (A_BITS_P - 1)
//
//       Stage B — c0 · Jx_i
//         c0jx_full [ACCUM_W_P+XY_FRAC_P] = Jx_i × zero_ext(c0_fp)
//         c0jx      [XY_W_P]              = c0jx_full >>> XY_FRAC_P
//
//       Stage C — Δy = dt · (−a·x + c0·Jx)
//         neg_ax    [XY_W_P]              = −ax_shifted
//         force     [XY_W_P]              = neg_ax + c0jx
//         dy_full   [2*XY_W_P]            = force × zero_ext(dt_fp)
//         delta_y   [XY_W_P]              = dy_full >>> XY_FRAC_P
//
//       Stage D — integrate y
//         y_pre     [XY_W_P]              = y_i + delta_y
//
//       Stage E — integrate x
//         dx_full   [2*XY_W_P]            = y_pre × zero_ext(dt_fp)
//         x_pre     [XY_W_P]              = x_i + (dx_full >>> XY_FRAC_P)
//
//       Stage F — inelastic wall
//         wall_hit_pos                    = x_pre >  ONE_FP
//         wall_hit_neg                    = x_pre < −ONE_FP
//         x_clamped [XY_W_P]             = mux(wall, ±ONE_FP, x_pre)
//         y_clamped [XY_W_P]             = mux(wall, 0,        y_pre)
//
//     The always_ff samples x_clamped → x_i_new, y_clamped → y_i_new
//     on posedge clk when start=1.
// =============================================================================
module update_unit
  import dsb_pkg::*;
#(
  parameter int unsigned XY_W_P    = XY_W,
  parameter int unsigned XY_FRAC_P = XY_FRAC,
  parameter int unsigned A_BITS_P  = A_BITS,
  parameter int unsigned ACCUM_W_P = ACCUM_W,
  parameter int unsigned PROD_W_P  = PROD_W
)(
  input  logic                          clk,
  input  logic                          start,
  input  logic signed [XY_W_P-1:0]     x_i,
  input  logic signed [XY_W_P-1:0]     y_i,
  input  logic signed [ACCUM_W_P-1:0]  Jx_i,
  input  logic        [A_BITS_P-1:0]   a_m,
  input  logic        [XY_FRAC_P-1:0]  dt_fp,
  input  logic        [XY_FRAC_P-1:0]  c0_fp,
  output logic signed [XY_W_P-1:0]     x_i_new,
  output logic signed [XY_W_P-1:0]     y_i_new,
  output logic                          done
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
  // Jx_i is a PLAIN SIGNED INTEGER count of J-unit contributions: it is the
  // direct accumulation of IC_BITS 2's-complement values from dotprod_phase2
  // (e.g. four neighbours all J=+1 gives Jx_i = +4, not scaled by 2^XY_FRAC).
  //
  // c0_fp is c0_real encoded as Q0.XY_FRAC_P, i.e. c0_fp = c0_real × 2^XY_FRAC.
  //
  // The fixed-point product we need is:
  //   c0_real × Jx_i   expressed as Q(XY_INT).(XY_FRAC) fixed-point
  //   = (c0_fp / 2^XY_FRAC) × Jx_i × 2^XY_FRAC
  //   = Jx_i × c0_fp                    ← NO division / shift needed
  //
  // This product can be large: for N=8, max|J|=127, Jx_max=1016, c0_fp=8192
  //   → c0jx_full_max = 1016 × 8192 = 8,323,072 >> XY_W signed max (32767).
  // Therefore c0jx_full MUST stay in its wide ACCUM_W_P+XY_FRAC_P register;
  // it is NOT truncated to XY_W_P here.  The truncation happens only after
  // the dt multiply in stage C below.
  logic signed [ACCUM_W_P+XY_FRAC_P-1:0] c0jx_full;
 
  always_comb begin
    c0jx_full = Jx_i * $signed({1'b0, c0_fp});   // wide — no shift
  end
 
  // ── Stage C: Δy = dt · (−a·x + c0·Jx) ───────────────────────────────────
  // Force is formed in the WIDE domain (ACCUM_W_P+XY_FRAC_P bits) to prevent
  // the c0·Jx term from overflowing before the dt multiply.
  //
  //   force_wide = −ax_shifted + c0jx_full
  //                (ax_shifted is XY_W_P wide — sign-extend before adding)
  //
  //   dy_full    = force_wide × dt_fp
  //                (dt_fp is Q0.XY_FRAC_P, so the product has an extra
  //                 factor of 2^XY_FRAC that we remove with one right-shift)
  //
  //   delta_y    = dy_full >>> XY_FRAC_P,  truncated to XY_W_P
  //
  // Width accounting (with default parameters):
  //   force_wide : ACCUM_W + XY_FRAC = 38 bits (signed)
  //   dy_full    : 38 + XY_FRAC = 52 bits (signed)
  //   delta_y    : XY_W bits after >>XY_FRAC
  localparam int unsigned FORCE_W = ACCUM_W_P + XY_FRAC_P;      // 38
  localparam int unsigned DY_W    = ACCUM_W_P + 2*XY_FRAC_P;    // 52
 
  logic signed [XY_W_P-1:0]   neg_ax;
  logic signed [FORCE_W-1:0]  force_wide;
  logic signed [DY_W-1:0]     dy_full;
  logic signed [XY_W_P-1:0]   delta_y;
 
  always_comb begin
    neg_ax     = -ax_shifted;
    // Sign-extend neg_ax to FORCE_W before adding the wide c0jx_full
    force_wide = $signed({{(FORCE_W-XY_W_P){neg_ax[XY_W_P-1]}}, neg_ax})
                 + $signed(c0jx_full);
    dy_full    = $signed(force_wide) * $signed({1'b0, dt_fp});
    delta_y    = XY_W_P'(dy_full >>> XY_FRAC_P);
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
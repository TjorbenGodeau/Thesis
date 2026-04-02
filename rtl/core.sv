module core
    import dsb_pkg::*;
    #(
    parameter int unsigned N_P       = N,
    parameter int unsigned IC_BITS_P = IC_BITS,
    parameter int unsigned XY_W_P    = XY_W,
    parameter int unsigned XY_FRAC_P = XY_FRAC,
    parameter int unsigned A_BITS_P  = A_BITS,
    parameter int unsigned ACCUM_W_P = ACCUM_W,
    parameter int unsigned PROD_W_P  = PROD_W
    )(
    input  logic                        clk,
    input  logic                        rst_n,
    input  logic                        start,
    // J row + neighbour sign write port
    input  logic                        wwl,
    input  logic [N_P*IC_BITS_P-1:0]   wbl_J,
    input  logic [N_P-1:0]              wbl_signs,
    // Current oscillator state
    input  logic signed [XY_W_P-1:0]   x_i,
    input  logic signed [XY_W_P-1:0]   y_i,
    // Schedule and constants
    input  logic [A_BITS_P-1:0]         a_m,
    input  logic [XY_FRAC_P-1:0]        dt_fp,
    input  logic [XY_FRAC_P-1:0]        c0_fp,
    // Outputs
    output logic signed [XY_W_P-1:0]   x_i_new,
    output logic signed [XY_W_P-1:0]   y_i_new,
    output logic                         sign_xi_new,  // sign(x_i_new)
    output logic                         done
    );

    // sign convention: 1 = positive, 0 = negative
    logic sign_xi_cur;
    assign sign_xi_cur = ~x_i[XY_W_P-1];
    assign sign_xi_new = ~x_i_new[XY_W_P-1];

    // --- tile interface ---
    dsb_tile_if #(
        .N_P (N_P),
        .IC_BITS_P (IC_BITS_P)
    ) tile_bus ();

    logic precharge_r;

    compute_tile #(
        .N_P (N_P),
        .IC_BITS_P (IC_BITS_P)
    ) u_tile (
        .clk (clk),
        .precharge (precharge_r),
        .wwl (wwl),
        .wbl_J (wbl_J),
        .wbl_signs (wbl_signs),
        .sign_xi (sign_xi_cur),
        .tile_bus (tile_bus.tile_out)
    );

    // --- Phase 1 ---
    logic capture_r;
    logic [N_P*IC_BITS_P-1:0] xnor_J;
    logic [N_P-1:0] sign_eq;
    
    xnor_phase1 #(
        .N_P (N_P), 
        .IC_BITS_P (IC_BITS_P)
    ) u_ph1 (
        .clk (clk),
        .capture (capture_r),
        .tile_bus (tile_bus.ph1_in),
        .xnor_J (xnor_J),
        .sign_eq (sign_eq)
    );

    // --- Phase 2 ---
    logic ph2_start;
    logic signed [ACCUM_W_P-1:0] Jx_i;
    logic ph2_done;
    logic sign_xi_r;  // latched at start
    
    dotprod_phase2 #(
        .N_P(N_P), 
        .IC_BITS_P(IC_BITS_P), 
        .ACCUM_W_P(ACCUM_W_P)
    ) u_ph2 (
        .clk (clk),
        .start (ph2_start),
        .sign_xi (sign_xi_r),
        .xnor_J (xnor_J),
        .sign_eq (sign_eq),
        .Jx_i (Jx_i),
        .done (ph2_done)
    );

    // --- Phase 3 ---
    logic ph3_start;
    logic signed [XY_W_P-1:0] x_new_w, y_new_w;
    logic ph3_done;
    logic signed [XY_W_P-1:0] x_i_r, y_i_r;  // latched at start
    
    update_unit #(
        .XY_W_P(XY_W_P), 
        .XY_FRAC_P(XY_FRAC_P), 
        .A_BITS_P(A_BITS_P),
        .ACCUM_W_P(ACCUM_W_P), 
        .PROD_W_P(PROD_W_P)
    ) u_ph3 (
        .clk (clk),
        .start (ph3_start),
        .x_i (x_i_r),
        .y_i (y_i_r),
        .Jx_i (Jx_i),
        .a_m (a_m),
        .dt_fp (dt_fp),
        .c0_fp (c0_fp),
        .x_i_new (x_new_w),
        .y_i_new (y_new_w),
        .done (ph3_done)
    );

    // --- FSM ---
    core_state_t state;
    
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state <= S_IDLE;
            precharge_r <= 1'b0;
            capture_r <= 1'b0;
            ph2_start <= 1'b0;
            ph3_start <= 1'b0;
            sign_xi_r <= 1'b0;
            x_i_r <= '0;
            y_i_r <= '0;
            x_i_new <= '0;
            y_i_new <= '0;
            done <= 1'b0;
        end else begin
            precharge_r <= 1'b0;
            capture_r <= 1'b0;
            ph2_start <= 1'b0;
            ph3_start <= 1'b0;
            done <= 1'b0;
    
            unique case (state)
                S_IDLE: begin
                    if (start) begin
                        x_i_r <= x_i;
                        y_i_r <= y_i;
                        sign_xi_r <= sign_xi_cur;
                        state <= S_PRECHARGE;
                    end
                end
                S_PRECHARGE: begin
                    precharge_r <= 1'b1;
                    state <= S_PHASE1;
                end
                S_PHASE1: begin
                    // RWL = sign_xi_r drives all bitcells; capture RBL output
                    capture_r <= 1'b1;
                    state <= S_PHASE2;
                end
                S_PHASE2: begin
                    ph2_start <= 1'b1;
                    state <= S_PHASE3;
                end
                S_PHASE3: begin
                    if (ph2_done) begin
                        ph3_start <= 1'b1;
                        state <= S_IDLE;
                    end
                end
                default: state <= S_IDLE;
            endcase
    
            // Latch phase-3 result whenever it finishes
            if (ph3_done) begin
                x_i_new <= x_new_w;
                y_i_new <= y_new_w;
                done <= 1'b1;
            end
        end
    end
endmodule
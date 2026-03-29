module dsb_array
    import dsb_pkg::*;
    #(
        parameter int unsigned N_P = N,
        parameter int unsigned IC_BITS_P = IC_BITS,
        parameter int unsigned XY_W_P = XY_W,
        parameter int unsigned XY_FRAC_P = XY_FRAC,
        parameter int unsigned A_BITS_P = A_BITS,
        parameter int unsigned STEP_W_P = STEP_W,
        parameter int unsigned ACCUM_W_P = ACCUM_W,
        parameter int unsigned PROD_W_P = PROD_W
    )(
        input  logic clk,
        input  logic rst_n,
        // ── Load interface ──────────────────────────────────────────────────────
        input  logic [$clog2(N_P)-1:0] wr_row,
        input  logic wr_en,
        input  logic [N_P*IC_BITS_P-1:0] wr_J_row,
        input  logic [$clog2(N_P)-1:0] wr_xy_idx,
        input  logic wr_xy_en,
        input  logic signed [XY_W_P-1:0] wr_x,
        input  logic signed [XY_W_P-1:0] wr_y,
        // ── Run control ─────────────────────────────────────────────────────────
        input  logic run,
        input  logic [STEP_W_P-1:0] Nstep,
        input  logic [A_BITS_P-1:0] a0_fp,
        input  logic [XY_FRAC_P-1:0] dt_fp,
        input  logic [XY_FRAC_P-1:0] c0_fp,
        // ── Outputs ─────────────────────────────────────────────────────────────
        output logic [N_P*XY_W_P-1:0] x_out,
        output logic [N_P*XY_W_P-1:0] y_out,
        output logic [N_P-1:0] signs_out,
        output logic step_done,
        output logic schedule_done
    );
 
    // ── State registers ────────────────────────────────────────────────────────
    logic signed [XY_W_P-1:0] x_reg [N_P];
    logic signed [XY_W_P-1:0] y_reg [N_P];
    logic [N_P*IC_BITS_P-1:0] J_table [N_P];
    
    // Unpack outputs
    for (genvar gi = 0; gi < N_P; gi++) begin : gen_out
        assign x_out[gi*XY_W_P +: XY_W_P] = x_reg[gi];
        assign y_out[gi*XY_W_P +: XY_W_P] = y_reg[gi];
        assign signs_out[gi]               = ~x_reg[gi][XY_W_P-1];
    end
    
    // ── Schedule ───────────────────────────────────────────────────────────────
    logic [STEP_W_P-1:0] m;
    logic [A_BITS_P-1:0] a_m;
    
    schedule #(
        .A_BITS_P (A_BITS_P), 
        .STEP_W_P (STEP_W_P)
    ) u_sched (
        .m (m),
        .Nstep (Nstep),
        .a0_fp (a0_fp),
        .a_m (a_m)
    );
    
    // ── Core ───────────────────────────────────────────────────────────────────
    logic [$clog2(N_P)-1:0] cur_osc;
    logic core_start;
    logic core_wwl;
    logic [N_P*IC_BITS_P-1:0] core_wbl_J;
    logic [N_P-1:0] core_wbl_signs;
    logic signed [XY_W_P-1:0] core_x_new, core_y_new;
    logic core_sign_new;
    logic core_done;
    
    core #(
        .N_P (N_P), 
        .IC_BITS_P (IC_BITS_P), 
        .XY_W_P (XY_W_P),
        .XY_FRAC_P (XY_FRAC_P), 
        .A_BITS_P (A_BITS_P),
        .ACCUM_W_P (ACCUM_W_P), 
        .PROD_W_P (PROD_W_P)
    ) u_core (
        .clk (clk),
        .rst_n (rst_n),
        .start (core_start),
        .wwl (core_wwl),
        .wbl_J (core_wbl_J),
        .wbl_signs (core_wbl_signs),
        .x_i (x_reg[cur_osc]),
        .y_i (y_reg[cur_osc]),
        .a_m (a_m),
        .dt_fp (dt_fp),
        .c0_fp (c0_fp),
        .x_i_new (core_x_new),
        .y_i_new (core_y_new),
        .sign_xi_new (core_sign_new),
        .done (core_done)
    );
    
    // ── Controller FSM ─────────────────────────────────────────────────────────
    ctrl_state_t ctrl;
    logic running;
    
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
        running <= 1'b0;
        m <= '0;
        cur_osc <= '0;
        step_done <= 1'b0;
        schedule_done <= 1'b0;
        core_start <= 1'b0;
        core_wwl <= 1'b0;
        ctrl <= CTRL_IDLE;
        for (int i = 0; i < N_P; i++) begin
            // Small non-zero initial x ≈ (i+1)×2^(−4) in fixed-point
            x_reg[i] <= XY_W_P'((i + 1) << (XY_FRAC_P - 4));
            y_reg[i] <= '0;
        end
        end else begin
            step_done <= 1'b0;
            core_start <= 1'b0;
            core_wwl <= 1'b0;
        
            // ── Load J rows ──────────────────────────────────────────────────────
            if (wr_en && !running) J_table[wr_row] <= wr_J_row;
        
            // ── Load initial x, y ─────────────────────────────────────────────────
            if (wr_xy_en && !running) begin
                x_reg[wr_xy_idx] <= wr_x;
                y_reg[wr_xy_idx] <= wr_y;
            end
        
            // ── Start schedule ────────────────────────────────────────────────────
            if (run && !running && !schedule_done) begin
                running <= 1'b1;
                m <= '0;
                ctrl <= CTRL_LOAD;
            end
        
            // ── Running ───────────────────────────────────────────────────────────
            if (running) begin
                unique case (ctrl)
                    CTRL_LOAD: begin
                        // Push J row and current sign snapshot into tile
                        core_wwl <= 1'b1;
                        core_wbl_J <= J_table[cur_osc];
                        core_wbl_signs <= signs_out;
                        ctrl <= CTRL_RUN;
                    end
                    CTRL_RUN: begin
                        core_start <= 1'b1;
                        ctrl <= CTRL_IDLE;
                    end
                    CTRL_IDLE: begin
                        if (core_done) begin
                            x_reg[cur_osc] <= core_x_new;
                            y_reg[cur_osc] <= core_y_new;
                
                            if (cur_osc == N_P - 1) begin
                                // All oscillators updated → dSB step complete
                                cur_osc <= '0;
                                step_done <= 1'b1;
                                if (m == Nstep - 1) begin
                                schedule_done <= 1'b1;
                                running <= 1'b0;
                                end else begin
                                m <= m + 1;
                                ctrl <= CTRL_LOAD;
                                end
                            end else begin
                                cur_osc <= cur_osc + 1;
                                ctrl <= CTRL_LOAD;
                            end
                        end
                    end
                    default: ctrl <= CTRL_IDLE;
                endcase
            end
        end
    end
 
endmodule
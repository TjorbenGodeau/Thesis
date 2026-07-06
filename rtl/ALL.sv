// =============================================================================
// sparse_dsb_top.sv
//
// Single-file synthesis-friendly version of the sparse dSB Ising machine.
// Combines all modules and packages into one self-contained file.
//
// Generated from:
//   allModules.sv, dsb_sparse_pkg.sv, csr_index_store.sv,
//   dotprod_phase2_sparse.sv, main_memory.sv, sparse_compute_memory.sv,
//   sparse_core.sv, sparse_dsb_array.sv, sparse_load_controller.sv,
//   sparsity_scanner.sv
// =============================================================================

// =============================================================================
// Package: dsb_pkg
// =============================================================================
package dsb_pkg;
    // --- Problem size ---
    parameter int unsigned N = 20000;   // max oscillators -- sized for the largest
                                         // Gset instance (G77/G81 = 20,000 nodes).
                                         // Smaller instances use the first
                                         // active_n of these (see sparse_dsb_array).
    parameter int unsigned IC_BITS = 7; // J coefficient bit width

    // --- Fixed-point format for x and y ---
    parameter int unsigned XY_INT = 2;
    parameter int unsigned XY_FRAC = 14;
    parameter int unsigned XY_W = XY_INT + XY_FRAC;

    // --- Schedule coefficient a(m) ---
    parameter int unsigned A_BITS = 16;

    // --- Step counter ---
    parameter int unsigned STEP_W = 12;     // supports Nstep up to 2^STEP_W - 1

    // --- JX accumulator (must hold N * max|J|) ---
    parameter int unsigned ACCUM_W = 24;

    // --- a*x product width ---
    parameter int unsigned PROD_W = XY_W + A_BITS;

    // --- Fixed-point unity for x ---
    parameter logic signed [XY_W-1:0] ONE_FP = XY_W'(1 << XY_FRAC);

    // --- Core FSM states ---
    typedef enum logic [2:0] {
        S_IDLE = 3'd0,
        S_PRECHARGE = 3'd1,
        S_PHASE1 = 3'd2,
        S_PHASE2 = 3'd3,
        S_PHASE3 = 3'd4
    } core_state_t;

    // --- Array controller FSM states ---
    typedef enum logic [1:0] {
        CTRL_IDLE = 2'd0,
        CTRL_LOAD = 2'd1,
        CTRL_RUN = 2'd2
    } ctrl_state_t;
endpackage

// =============================================================================
// Package: dsb_sparse_pkg
// =============================================================================
package dsb_sparse_pkg;
    import dsb_pkg::*;

    // -------------------------------------------------------------------------
    // Sparse compute memory geometry
    // -------------------------------------------------------------------------
    parameter int unsigned K_MAX = 32;

    // -------------------------------------------------------------------------
    // Edge-list capacity
    // -------------------------------------------------------------------------
    parameter int unsigned MAX_EDGES = 50_000;
    parameter int unsigned MAX_NNZ   = 2 * MAX_EDGES;

    // -------------------------------------------------------------------------
    // Derived widths
    // -------------------------------------------------------------------------
    parameter int unsigned COL_IDX_W = $clog2(N);

    // -------------------------------------------------------------------------
    // Main memory layout — EDGE LIST based
    // -------------------------------------------------------------------------
    parameter int unsigned EDGE_ROW_W   = COL_IDX_W;
    parameter int unsigned EDGE_COL_W   = COL_IDX_W;
    parameter int unsigned EDGE_WT_W    = IC_BITS;
    parameter int unsigned EDGE_WORD_W  = EDGE_ROW_W + EDGE_COL_W + EDGE_WT_W;
    parameter int unsigned MAIN_WORD_W  = (EDGE_WORD_W > XY_W) ? EDGE_WORD_W : XY_W;
    parameter int unsigned MAIN_MEM_WORDS = MAX_EDGES + 2 * N;
    parameter int unsigned MAIN_ADDR_W    = $clog2(MAIN_MEM_WORDS);

    parameter int unsigned MAIN_EDGE_BASE = 0;
    parameter int unsigned MAIN_X_BASE    = MAX_EDGES;
    parameter int unsigned MAIN_Y_BASE    = MAX_EDGES + N;

    // -------------------------------------------------------------------------
    // CSR entry store address width
    // -------------------------------------------------------------------------
    parameter int unsigned CSR_ADDR_W  = $clog2(MAX_NNZ + 1);
    parameter int unsigned CSR_ENTRY_W = COL_IDX_W + MAIN_ADDR_W;

    parameter int unsigned ROW_COUNT_W = $clog2(MAX_NNZ + 1);
    parameter int unsigned ROW_PTR_W   = CSR_ADDR_W + ROW_COUNT_W;

    // -------------------------------------------------------------------------
    // Chunk counter
    // -------------------------------------------------------------------------
    parameter int unsigned MAX_ROW_NNZ = MAX_NNZ;
    parameter int unsigned MAX_CHUNKS  = (MAX_ROW_NNZ + K_MAX - 1) / K_MAX;
    parameter int unsigned CHUNK_CNT_W = $clog2(MAX_CHUNKS + 1);

    // -------------------------------------------------------------------------
    // FSM state types
    // -------------------------------------------------------------------------
    typedef enum logic [2:0] {
        SCAN_IDLE     = 3'd0,
        SCAN_READ     = 3'd1,
        SCAN_WAIT     = 3'd2,
        SCAN_WR_FWD   = 3'd3,
        SCAN_WR_REV   = 3'd4,
        SCAN_NEXTEDGE = 3'd5,
        SCAN_FINALIZE = 3'd6,
        SCAN_DONE     = 3'd7
    } scan_state_t;

    typedef enum logic [3:0] {
        SLC_IDLE       = 4'd0,
        SLC_LOAD_PTR   = 4'd1,
        SLC_WAIT_PTR   = 4'd2,
        SLC_CHUNK_START= 4'd3,
        SLC_LOAD_ENTRY = 4'd4,
        SLC_WAIT_ENTRY = 4'd5,
        SLC_READ_J     = 4'd6,
        SLC_WAIT_J     = 4'd7,
        SLC_WRITE_TILE = 4'd8,
        SLC_NEXT_ENTRY = 4'd9,
        SLC_TRIGGER_PH2= 4'd10,
        SLC_WAIT_PH2   = 4'd11,
        SLC_NEXT_CHUNK = 4'd12,
        SLC_DONE       = 4'd13
    } slc_state_t;
endpackage

// =============================================================================
// Module: dsb_8t_bitcell
// =============================================================================
module dsb_8t_bitcell(
    input logic clk,
    input logic wwl,
    input logic wbl,
    input logic rwl,
    input logic precharge,
    output logic rbl_out
);
    logic q, qb;

    always_ff @(posedge clk) begin
        if (wwl) begin
            q <= wbl;
            qb <= ~wbl;
        end
    end

    always_comb begin
        if (precharge) rbl_out = 1'b1;
        else if (wwl) rbl_out = 1'b1;
        else rbl_out = ~(q ^ rwl);
    end
endmodule

// =============================================================================
// Module: schedule
// =============================================================================
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

// =============================================================================
// Module: update_unit
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
    localparam int unsigned FORCE_W = ACCUM_W_P + XY_FRAC_P;
    localparam int unsigned DY_W    = ACCUM_W_P + 2*XY_FRAC_P;

    logic signed [PROD_W_P-1:0]   ax_full;
    logic signed [XY_W_P-1:0]     ax_shifted;
    logic signed [ACCUM_W_P+XY_FRAC_P-1:0] c0jx_full;
    logic signed [XY_W_P-1:0]   neg_ax;
    logic signed [FORCE_W-1:0]  force_wide;
    logic signed [DY_W-1:0]     dy_full;
    logic signed [XY_W_P-1:0]   delta_y;
    logic signed [XY_W_P-1:0]   y_pre;
    logic signed [2*XY_W_P-1:0] dx_full;
    logic signed [XY_W_P-1:0]   x_pre;
    logic wall_hit_pos;
    logic wall_hit_neg;
    logic wall_hit;
    logic signed [XY_W_P-1:0] x_clamped;
    logic signed [XY_W_P-1:0] y_clamped;

    always_comb begin
        ax_full   = $signed({{(PROD_W_P-XY_W_P){x_i[XY_W_P-1]}}, x_i})
                    * $signed({1'b0, a_m});
        ax_shifted = XY_W_P'(ax_full >>> (A_BITS_P - 1));
    end

    always_comb begin
        c0jx_full = Jx_i * $signed({1'b0, c0_fp});
    end

    always_comb begin
        neg_ax     = -ax_shifted;
        force_wide = $signed({{(FORCE_W-XY_W_P){neg_ax[XY_W_P-1]}}, neg_ax})
                     + $signed(c0jx_full);
        dy_full    = $signed(force_wide) * $signed({1'b0, dt_fp});
        delta_y    = XY_W_P'(dy_full >>> XY_FRAC_P);
    end

    always_comb begin
        y_pre = y_i + delta_y;
    end

    always_comb begin
        dx_full = $signed(y_pre) * $signed({1'b0, dt_fp});
        x_pre   = x_i + XY_W_P'(dx_full >>> XY_FRAC_P);
    end

    always_comb begin
        wall_hit_pos = (x_pre >  $signed(ONE_FP));
        wall_hit_neg = (x_pre < -$signed(ONE_FP));
        wall_hit     = wall_hit_pos | wall_hit_neg;

        unique case ({wall_hit_pos, wall_hit_neg})
            2'b10   : x_clamped =  ONE_FP;
            2'b01   : x_clamped = -ONE_FP;
            default : x_clamped =  x_pre;
        endcase

        y_clamped = wall_hit ? '0 : y_pre;
    end

    always_ff @(posedge clk) begin
        done <= 1'b0;
        if (start) begin
            x_i_new <= x_clamped;
            y_i_new <= y_clamped;
            done    <= 1'b1;
        end
    end
endmodule

// =============================================================================
// Module: main_memory
// =============================================================================
module main_memory
    import dsb_pkg::*;
    import dsb_sparse_pkg::*;
    #(
        parameter int unsigned WORDS_P    = MAIN_MEM_WORDS,
        parameter int unsigned WORD_W_P   = MAIN_WORD_W,
        parameter int unsigned ADDR_W_P   = MAIN_ADDR_W
    )(
        input  logic                  clk,
        input  logic                  rst_n,
        input  logic                  wr_en,
        input  logic [ADDR_W_P-1:0]   wr_addr,
        input  logic [WORD_W_P-1:0]   wr_data,
        input  logic                  rd_en,
        input  logic [ADDR_W_P-1:0]   rd_addr,
        output logic [WORD_W_P-1:0]   rd_data,
        output logic                  rd_valid
    );
    logic [WORD_W_P-1:0] mem [0:WORDS_P-1];

    always_ff @(posedge clk) begin
        if (wr_en)
            mem[wr_addr] <= wr_data;
    end

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            rd_data  <= '0;
            rd_valid <= 1'b0;
        end else begin
            rd_valid <= 1'b0;
            if (rd_en) begin
                rd_data  <= mem[rd_addr];
                rd_valid <= 1'b1;
            end
        end
    end
endmodule

// =============================================================================
// Module: csr_index_store
// =============================================================================
module csr_index_store
    import dsb_pkg::*;
    import dsb_sparse_pkg::*;
    #(
        parameter int unsigned N_P          = N,
        parameter int unsigned MAX_NNZ_P    = MAX_NNZ,
        parameter int unsigned COL_IDX_W_P  = COL_IDX_W,
        parameter int unsigned MAIN_ADDR_WP = MAIN_ADDR_W,
        parameter int unsigned CSR_ADDR_W_P = CSR_ADDR_W,
        parameter int unsigned CSR_ENTRY_WP = CSR_ENTRY_W,
        parameter int unsigned ROW_PTR_W_P  = ROW_PTR_W,
        parameter int unsigned ROW_COUNT_WP = ROW_COUNT_W
    )(
        input  logic clk,
        input  logic rst_n,
        input  logic                    rp_wr_en,
        input  logic [COL_IDX_W_P-1:0]  rp_wr_row,
        input  logic [ROW_PTR_W_P-1:0]  rp_wr_data,
        input  logic                    rp_rd_en,
        input  logic [COL_IDX_W_P-1:0]  rp_rd_row,
        output logic [ROW_PTR_W_P-1:0]  rp_rd_data,
        output logic                    rp_rd_valid,
        input  logic                    en_wr_en,
        input  logic [CSR_ADDR_W_P-1:0] en_wr_addr,
        input  logic [CSR_ENTRY_WP-1:0] en_wr_data,
        input  logic                    en_rd_en,
        input  logic [CSR_ADDR_W_P-1:0] en_rd_addr,
        output logic [CSR_ENTRY_WP-1:0] en_rd_data,
        output logic                    en_rd_valid
    );
    logic [ROW_PTR_W_P-1:0] row_ptr_mem [0:N_P-1];
    logic [CSR_ENTRY_WP-1:0] entry_mem [0:MAX_NNZ_P-1];

    always_ff @(posedge clk) begin
        if (rp_wr_en)
            row_ptr_mem[rp_wr_row] <= rp_wr_data;
    end

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            rp_rd_data  <= '0;
            rp_rd_valid <= 1'b0;
        end else begin
            rp_rd_valid <= 1'b0;
            if (rp_rd_en) begin
                rp_rd_data  <= row_ptr_mem[rp_rd_row];
                rp_rd_valid <= 1'b1;
            end
        end
    end

    always_ff @(posedge clk) begin
        if (en_wr_en)
            entry_mem[en_wr_addr] <= en_wr_data;
    end

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            en_rd_data  <= '0;
            en_rd_valid <= 1'b0;
        end else begin
            en_rd_valid <= 1'b0;
            if (en_rd_en) begin
                en_rd_data  <= entry_mem[en_rd_addr];
                en_rd_valid <= 1'b1;
            end
        end
    end
endmodule

// =============================================================================
// Module: dotprod_phase2_sparse
// =============================================================================
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
            automatic logic signed [ACCUM_W_P-1:0] chunk_sum;
            chunk_sum = '0;

            for (int k = 0; k < K_MAX_P; k++) begin
                if (k < int'(valid_count)) begin
                    automatic logic [IC_BITS_P-1:0] xnor_word;
                    automatic logic [IC_BITS_P-1:0] corrected;
                    xnor_word = xnor_J[k*IC_BITS_P +: IC_BITS_P];
                    corrected  = sign_eq[k] ? xnor_word : ~xnor_word;
                    chunk_sum += ACCUM_W_P'(signed'(corrected));
                end
            end

            begin
                automatic logic signed [ACCUM_W_P-1:0] base;
                base = (accumulate ? Jx_accum : '0) + chunk_sum;

                if (last_chunk && !sign_xi)
                    Jx_accum <= base + ACCUM_W_P'(valid_count);
                else
                    Jx_accum <= base;
            end

            done <= 1'b1;
        end
    end

    assign Jx_i = Jx_accum;
endmodule

// =============================================================================
// Module: sparsity_scanner
// =============================================================================
module sparsity_scanner
    import dsb_pkg::*;
    import dsb_sparse_pkg::*;
    #(
        parameter int unsigned N_P          = N,
        parameter int unsigned IC_BITS_P    = IC_BITS,
        parameter int unsigned MAIN_ADDR_WP = MAIN_ADDR_W,
        parameter int unsigned MAIN_WORD_WP = MAIN_WORD_W,
        parameter int unsigned COL_IDX_W_P  = COL_IDX_W,
        parameter int unsigned CSR_ADDR_W_P = CSR_ADDR_W,
        parameter int unsigned CSR_ENTRY_WP = CSR_ENTRY_W,
        parameter int unsigned ROW_PTR_W_P  = ROW_PTR_W,
        parameter int unsigned ROW_COUNT_WP = ROW_COUNT_W
    )(
        input  logic clk,
        input  logic rst_n,
        input  logic start,
        input  logic [CSR_ADDR_W_P-1:0]  num_edges,
        input  logic [COL_IDX_W_P-1:0]   active_n,
        output logic                     mm_rd_en,
        output logic [MAIN_ADDR_WP-1:0]  mm_rd_addr,
        input  logic [MAIN_WORD_WP-1:0]  mm_rd_data,
        input  logic                     mm_rd_valid,
        output logic                     rp_wr_en,
        output logic [COL_IDX_W_P-1:0]   rp_wr_row,
        output logic [ROW_PTR_W_P-1:0]   rp_wr_data,
        output logic                     en_wr_en,
        output logic [CSR_ADDR_W_P-1:0]  en_wr_addr,
        output logic [CSR_ENTRY_WP-1:0]  en_wr_data,
        output logic                     done,
        output logic [CSR_ADDR_W_P-1:0]  nnz_total
    );
    logic [ROW_COUNT_WP-1:0] degree [0:N_P-1];
    logic [CSR_ADDR_W_P-1:0] rstart [0:N_P-1];
    logic [ROW_COUNT_WP-1:0] fill   [0:N_P-1];

    typedef enum logic [3:0] {
        ST_IDLE,
        ST_CLEAR,
        ST_P1_READ, ST_P1_WAIT, ST_P1_ACC,
        ST_FIN_INIT,
        ST_FIN_STEP,
        ST_FIN_WRROWPTR,
        ST_P2_READ, ST_P2_WAIT, ST_P2_PLACE_FWD, ST_P2_PLACE_REV,
        ST_DONE
    } scan2_state_t;

    scan2_state_t state;
    logic [CSR_ADDR_W_P-1:0] edge_idx;
    logic [COL_IDX_W_P-1:0]  clear_idx;
    logic [CSR_ADDR_W_P-1:0] running_sum;
    logic [COL_IDX_W_P-1:0] e_row, e_col;
    logic [MAIN_ADDR_WP-1:0] e_addr;

    function automatic logic [COL_IDX_W_P-1:0] unpack_row(logic [MAIN_WORD_WP-1:0] w);
        return w[MAIN_WORD_WP-1 -: COL_IDX_W_P];
    endfunction
    function automatic logic [COL_IDX_W_P-1:0] unpack_col(logic [MAIN_WORD_WP-1:0] w);
        return w[IC_BITS_P +: COL_IDX_W_P];
    endfunction

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state        <= ST_IDLE;
            edge_idx     <= '0;
            clear_idx    <= '0;
            running_sum  <= '0;
            done         <= 1'b0;
            nnz_total    <= '0;
            mm_rd_en     <= 1'b0;
            mm_rd_addr   <= '0;
            rp_wr_en     <= 1'b0;
            rp_wr_row    <= '0;
            rp_wr_data   <= '0;
            en_wr_en     <= 1'b0;
            en_wr_addr   <= '0;
            en_wr_data   <= '0;
        end else begin
            mm_rd_en <= 1'b0;
            rp_wr_en <= 1'b0;
            en_wr_en <= 1'b0;

            unique case (state)
                ST_IDLE: begin
                    if (start) begin
                        clear_idx <= '0;
                        done      <= 1'b0;
                        state     <= ST_CLEAR;
                    end
                end

                ST_CLEAR: begin
                    degree[clear_idx] <= '0;
                    fill[clear_idx]   <= '0;
                    if (clear_idx == COL_IDX_W_P'(active_n - 1)) begin
                        edge_idx <= '0;
                        state    <= ST_P1_READ;
                    end else begin
                        clear_idx <= clear_idx + 1;
                    end
                end

                ST_P1_READ: begin
                    if (edge_idx == num_edges) begin
                        clear_idx   <= '0;
                        running_sum <= '0;
                        state       <= ST_FIN_INIT;
                    end else begin
                        mm_rd_en   <= 1'b1;
                        mm_rd_addr <= MAIN_ADDR_WP'(MAIN_EDGE_BASE + edge_idx);
                        state      <= ST_P1_WAIT;
                    end
                end

                ST_P1_WAIT: state <= ST_P1_ACC;

                ST_P1_ACC: begin
                    if (mm_rd_valid) begin
                        automatic logic [COL_IDX_W_P-1:0] r, c;
                        r = unpack_row(mm_rd_data);
                        c = unpack_col(mm_rd_data);
                        degree[r] <= degree[r] + 1;
                        degree[c] <= degree[c] + 1;
                    end
                    edge_idx <= edge_idx + 1;
                    state    <= ST_P1_READ;
                end

                ST_FIN_INIT: state <= ST_FIN_STEP;

                ST_FIN_STEP: begin
                    rstart[clear_idx] <= running_sum;
                    running_sum       <= running_sum + CSR_ADDR_W_P'(degree[clear_idx]);
                    state             <= ST_FIN_WRROWPTR;
                end

                ST_FIN_WRROWPTR: begin
                    rp_wr_en   <= 1'b1;
                    rp_wr_row  <= clear_idx;
                    rp_wr_data <= {rstart[clear_idx], ROW_COUNT_WP'(degree[clear_idx])};

                    if (clear_idx == COL_IDX_W_P'(active_n - 1)) begin
                        nnz_total <= running_sum;
                        edge_idx  <= '0;
                        state     <= ST_P2_READ;
                    end else begin
                        clear_idx <= clear_idx + 1;
                        state     <= ST_FIN_STEP;
                    end
                end

                ST_P2_READ: begin
                    if (edge_idx == num_edges) begin
                        state <= ST_DONE;
                    end else begin
                        mm_rd_en   <= 1'b1;
                        mm_rd_addr <= MAIN_ADDR_WP'(MAIN_EDGE_BASE + edge_idx);
                        state      <= ST_P2_WAIT;
                    end
                end

                ST_P2_WAIT: state <= ST_P2_PLACE_FWD;

                ST_P2_PLACE_FWD: begin
                    if (mm_rd_valid) begin
                        automatic logic [COL_IDX_W_P-1:0] row_v, col_v;
                        automatic logic [CSR_ADDR_W_P-1:0] dest_addr;
                        row_v = unpack_row(mm_rd_data);
                        col_v = unpack_col(mm_rd_data);

                        e_row  <= row_v;
                        e_col  <= col_v;
                        e_addr <= MAIN_ADDR_WP'(MAIN_EDGE_BASE + edge_idx);

                        dest_addr  = CSR_ADDR_W_P'(rstart[row_v] + fill[row_v]);
                        en_wr_en   <= 1'b1;
                        en_wr_addr <= dest_addr;
                        en_wr_data <= {col_v, MAIN_ADDR_WP'(MAIN_EDGE_BASE + edge_idx)};
                        fill[row_v] <= fill[row_v] + 1;
                    end
                    state <= ST_P2_PLACE_REV;
                end

                ST_P2_PLACE_REV: begin
                    en_wr_en   <= 1'b1;
                    en_wr_addr <= CSR_ADDR_W_P'(rstart[e_col] + fill[e_col]);
                    en_wr_data <= {e_row, e_addr};
                    fill[e_col] <= fill[e_col] + 1;

                    edge_idx <= edge_idx + 1;
                    state    <= ST_P2_READ;
                end

                ST_DONE: begin
                    done <= 1'b1;
                    if (start) begin
                        clear_idx <= '0;
                        done      <= 1'b0;
                        state     <= ST_CLEAR;
                    end
                end

                default: state <= ST_IDLE;
            endcase
        end
    end
endmodule

// =============================================================================
// Module: sparse_compute_memory
// =============================================================================
module sparse_compute_memory
    import dsb_pkg::*;
    import dsb_sparse_pkg::*;
    #(
        parameter int unsigned K_MAX_P   = K_MAX,
        parameter int unsigned IC_BITS_P = IC_BITS,
        parameter int unsigned COL_W_P   = COL_IDX_W
    )(
        input  logic clk,
        input  logic                      wr_en,
        input  logic [$clog2(K_MAX)-1:0]  wr_slot,
        input  logic [COL_W_P-1:0]        wr_col_idx,
        input  logic [IC_BITS_P-1:0]      wr_J,
        input  logic                      wr_sign_j,
        input  logic                      rd_col_en,
        input  logic [$clog2(K_MAX)-1:0]  rd_slot,
        output logic [COL_W_P-1:0]        rd_col_idx,
        output logic                      rd_col_valid,
        input  logic                      precharge,
        input  logic [K_MAX_P-1:0]        rwl,
        output logic [K_MAX_P*IC_BITS_P-1:0] rbl_J,
        output logic [K_MAX_P-1:0]        rbl_signs
    );
    logic [COL_W_P-1:0] col_idx_mem [0:K_MAX_P-1];

    always_ff @(posedge clk) begin
        if (wr_en)
            col_idx_mem[wr_slot] <= wr_col_idx;
    end

    always_ff @(posedge clk) begin
        rd_col_valid <= 1'b0;
        if (rd_col_en) begin
            rd_col_idx   <= col_idx_mem[rd_slot];
            rd_col_valid <= 1'b1;
        end
    end

    for (genvar k = 0; k < K_MAX_P; k++) begin : gen_slot
        logic wwl_k;
        assign wwl_k = wr_en && (wr_slot == $clog2(K_MAX)'(k));

        for (genvar b = 0; b < IC_BITS_P; b++) begin : gen_bit
            dsb_8t_bitcell u_J (
                .clk       (clk),
                .wwl       (wwl_k),
                .wbl       (wr_J[b]),
                .rwl       (rwl[k]),
                .precharge (precharge),
                .rbl_out   (rbl_J[k*IC_BITS_P + b])
            );
        end

        dsb_8t_bitcell u_sign (
            .clk       (clk),
            .wwl       (wwl_k),
            .wbl       (wr_sign_j),
            .rwl       (rwl[k]),
            .precharge (precharge),
            .rbl_out   (rbl_signs[k])
        );
    end
endmodule

// =============================================================================
// Module: sparse_load_controller
// =============================================================================
module sparse_load_controller
    import dsb_pkg::*;
    import dsb_sparse_pkg::*;
    #(
        parameter int unsigned N_P          = N,
        parameter int unsigned K_MAX_P      = K_MAX,
        parameter int unsigned IC_BITS_P    = IC_BITS,
        parameter int unsigned COL_IDX_W_P  = COL_IDX_W,
        parameter int unsigned MAIN_ADDR_WP = MAIN_ADDR_W,
        parameter int unsigned MAIN_WORD_WP = MAIN_WORD_W,
        parameter int unsigned CSR_ADDR_W_P = CSR_ADDR_W,
        parameter int unsigned CSR_ENTRY_WP = CSR_ENTRY_W,
        parameter int unsigned ROW_PTR_W_P  = ROW_PTR_W,
        parameter int unsigned ROW_COUNT_WP = ROW_COUNT_W,
        parameter int unsigned K_CNT_W_P    = $clog2(K_MAX + 1),
        parameter int unsigned XY_W_P       = XY_W
    )(
        input  logic clk,
        input  logic rst_n,
        input  logic                        start,
        input  logic [COL_IDX_W_P-1:0]      osc_idx,
        output logic                        done,
        input  logic [N_P-1:0]              osc_signs,
        output logic                        rp_rd_en,
        output logic [COL_IDX_W_P-1:0]      rp_rd_row,
        input  logic [ROW_PTR_W_P-1:0]      rp_rd_data,
        input  logic                        rp_rd_valid,
        output logic                        en_rd_en,
        output logic [CSR_ADDR_W_P-1:0]     en_rd_addr,
        input  logic [CSR_ENTRY_WP-1:0]     en_rd_data,
        input  logic                        en_rd_valid,
        output logic                        mm_rd_en,
        output logic [MAIN_ADDR_WP-1:0]     mm_rd_addr,
        input  logic [MAIN_WORD_WP-1:0]     mm_rd_data,
        input  logic                        mm_rd_valid,
        output logic                        scm_wr_en,
        output logic [$clog2(K_MAX)-1:0]    scm_wr_slot,
        output logic [COL_IDX_W_P-1:0]      scm_wr_col_idx,
        output logic [IC_BITS_P-1:0]        scm_wr_J,
        output logic                        scm_wr_sign_j,
        output logic                        scm_precharge,
        output logic [K_MAX_P-1:0]          scm_rwl,
        output logic                        ph2_start,
        output logic                        ph2_accumulate,
        output logic                        ph2_last_chunk,
        output logic [K_CNT_W_P-1:0]        ph2_valid_count,
        output logic                        ph2_sign_xi,
        input  logic                        ph2_done
    );
    slc_state_t state;

    logic [COL_IDX_W_P-1:0]   cur_osc;
    logic [CSR_ADDR_W_P-1:0]  csr_start;
    logic [ROW_COUNT_WP-1:0]  row_count;
    logic [ROW_COUNT_WP-1:0]  entries_remaining;
    logic [CSR_ADDR_W_P-1:0]  entry_ptr;
    logic [$clog2(K_MAX)-1:0] slot_ptr;
    logic [K_CNT_W_P-1:0]     chunk_valid;
    logic                      first_chunk;

    logic [COL_IDX_W_P-1:0]   latch_col_j;
    logic [MAIN_ADDR_WP-1:0]  latch_main_addr;
    logic [IC_BITS_P-1:0]     latch_J;

    logic sign_xi;
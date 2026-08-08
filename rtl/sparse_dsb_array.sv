// =============================================================================
// sparse_dsb_array.sv
//
// Top-level module for the SPARSE dSB Ising machine.
// Replaces the original dsb_array, keeping its external interface identical
// so it can be swapped in without changing any upstream testbench or system.
//
// What is different from the original dsb_array
// -----------------------------------------------
//  REMOVED : J_table[]     (dense N×N row store inside the array)
//  REMOVED : wr_en / wr_J_row  (dense J push interface)
//  REMOVED : core          (replaced by sparse_core)
//
//  ADDED   : main_memory        (holds J matrix + x/y init, word-addressable)
//  ADDED   : sparsity_scanner   (startup FSM, builds CSR store once)
//  ADDED   : csr_index_store    (row_ptr table + flat entry array)
//  ADDED   : sparse_core        (SLC + sparse_compute_memory + dotprod_phase2_sparse
//                                + update_unit)
//  ADDED   : startup FSM (BOOT phase): waits for scan_done before accepting run
//
// External interface changes
// --------------------------
//  REMOVED : wr_en, wr_J_row    (J is loaded via mm_wr_* below instead)
//  ADDED   : mm_wr_en, mm_wr_addr, mm_wr_data  (word write port into main_memory)
//  ADDED   : active_n   (oscillators actually used by this Gset instance)
//  ADDED   : num_edges  (edges actually loaded for this Gset instance)
//
// All other ports (wr_xy_en, wr_xy_idx, wr_x, wr_y, run, Nstep, a0_fp, dt_fp,
// c0_fp, x_out, y_out, signs_out, step_done, schedule_done) are IDENTICAL.
//
// Startup sequence (v3 — edge-list based)
// -----------------------------------------
//  1. Write the edge list:  for each edge e = 0..num_edges-1,
//       mm_wr_en=1, mm_wr_addr=MAIN_EDGE_BASE+e,
//       mm_wr_data = {row[e], col[e], weight[e]}   (packed per dsb_sparse_pkg)
//  2. Write x/y init:  wr_xy_en=1 (unchanged from original dsb_array)
//  3. Set active_n = number of nodes in this instance, num_edges = edge count
//  4. Assert run=1.  The internal boot FSM automatically fires the
//     edge-list-ingest scanner and waits for scan_done before starting
//     the dSB sweep over oscillators 0..active_n-1.
// =============================================================================

`timescale 1ns/1ps
module sparse_dsb_array
    import dsb_pkg::*;
    import dsb_sparse_pkg::*;
#(
    parameter int unsigned N_P        = N,
    parameter int unsigned IC_BITS_P  = IC_BITS,
    parameter int unsigned XY_W_P     = XY_W,
    parameter int unsigned XY_FRAC_P  = XY_FRAC,
    parameter int unsigned A_BITS_P   = A_BITS,
    parameter int unsigned STEP_W_P   = STEP_W,
    parameter int unsigned ACCUM_W_P  = ACCUM_W,
    parameter int unsigned PROD_W_P   = PROD_W
)(
    input  logic clk,
    input  logic rst_n,

    // ── Main memory write port (replaces wr_en / wr_J_row) ───────────────────
    input  logic                    mm_wr_en,
    input  logic [MAIN_ADDR_W-1:0]  mm_wr_addr,
    input  logic [MAIN_WORD_W-1:0]  mm_wr_data,

    // ── Problem-instance sizing (lets one build run any Gset instance) ───────
    // active_n   : number of oscillators actually used by THIS instance
    //              (<= N_P; e.g. 800 for G1, 20000 for G81)
    // num_edges  : number of edges present in the edge-list region of
    //              main_memory for THIS instance
    input  logic [$clog2(N_P)-1:0]  active_n,
    input  logic [CSR_ADDR_W-1:0]   num_edges,

    // ── Initial x/y write port (identical to original dsb_array) ─────────────
    input  logic [$clog2(N_P)-1:0]  wr_xy_idx,
    input  logic                    wr_xy_en,
    input  logic signed [XY_W_P-1:0] wr_x,
    input  logic signed [XY_W_P-1:0] wr_y,

    // ── Run control (identical to original dsb_array) ─────────────────────────
    input  logic                    run,
    input  logic [STEP_W_P-1:0]     Nstep,
    input  logic [A_BITS_P-1:0]     a0_fp,
    input  logic [XY_FRAC_P-1:0]    dt_fp,
    input  logic [XY_FRAC_P-1:0]    c0_fp,

    // ── Outputs (identical to original dsb_array) ─────────────────────────────
    output logic [N_P*XY_W_P-1:0]  x_out,
    output logic [N_P*XY_W_P-1:0]  y_out,
    output logic [N_P-1:0]          signs_out,
    output logic                    step_done,
    output logic                    schedule_done,

    // ── Extra: scan status ────────────────────────────────────────────────────
    output logic                    scan_done       // high once sparsity scan complete
);

    // =========================================================================
    // State registers  (same as original dsb_array)
    // =========================================================================
    logic signed [XY_W_P-1:0] x_reg [N_P];
    logic signed [XY_W_P-1:0] y_reg [N_P];

    // Unpack outputs (identical logic to original)
    for (genvar gi = 0; gi < N_P; gi++) begin : gen_out
        assign x_out[gi*XY_W_P +: XY_W_P] = x_reg[gi];
        assign y_out[gi*XY_W_P +: XY_W_P] = y_reg[gi];
        assign signs_out[gi]               = ~x_reg[gi][XY_W_P-1];
    end

    // =========================================================================
    // Schedule  (identical to original dsb_array)
    // =========================================================================
    logic [STEP_W_P-1:0] m;
    logic [A_BITS_P-1:0] a_m;

    schedule #(
        .A_BITS_P(A_BITS_P),
        .STEP_W_P(STEP_W_P)
    ) u_sched (
        .m      (m),
        .Nstep  (Nstep),
        .a0_fp  (a0_fp),
        .a_m    (a_m)
    );

    // =========================================================================
    // main_memory
    // =========================================================================
    logic                   mm_rd_en;
    logic [MAIN_ADDR_W-1:0] mm_rd_addr;
    logic [MAIN_WORD_W-1:0] mm_rd_data;
    logic                   mm_rd_valid;

    // Read bus mux: scanner owns it during BOOT, sparse_core owns it during RUN
    logic                   sc_rd_en;
    logic [MAIN_ADDR_W-1:0] sc_rd_addr;
    logic                   core_rd_en;
    logic [MAIN_ADDR_W-1:0] core_rd_addr;

    assign mm_rd_en   = scan_done ? core_rd_en   : sc_rd_en;
    assign mm_rd_addr = scan_done ? core_rd_addr  : sc_rd_addr;

    main_memory u_mm (
        .clk      (clk),
        .rst_n    (rst_n),
        .wr_en    (mm_wr_en),
        .wr_addr  (mm_wr_addr),
        .wr_data  (mm_wr_data),
        .rd_en    (mm_rd_en),
        .rd_addr  (mm_rd_addr),
        .rd_data  (mm_rd_data),
        .rd_valid (mm_rd_valid)
    );

    // =========================================================================
    // csr_index_store
    // =========================================================================
    // Write ports driven by scanner
    logic                   rp_wr_en;
    logic [COL_IDX_W-1:0]   rp_wr_row;
    logic [ROW_PTR_W-1:0]   rp_wr_data;
    logic                   en_wr_en;
    logic [CSR_ADDR_W-1:0]  en_wr_addr;
    logic [CSR_ENTRY_W-1:0] en_wr_data;
    // Read ports driven by sparse_core
    logic                   rp_rd_en;
    logic [COL_IDX_W-1:0]   rp_rd_row;
    logic [ROW_PTR_W-1:0]   rp_rd_data;
    logic                   rp_rd_valid;
    logic                   en_rd_en;
    logic [CSR_ADDR_W-1:0]  en_rd_addr;
    logic [CSR_ENTRY_W-1:0] en_rd_data;
    logic                   en_rd_valid;

    csr_index_store u_csr (
        .clk         (clk),
        .rst_n       (rst_n),
        .rp_wr_en    (rp_wr_en),
        .rp_wr_row   (rp_wr_row),
        .rp_wr_data  (rp_wr_data),
        .rp_rd_en    (rp_rd_en),
        .rp_rd_row   (rp_rd_row),
        .rp_rd_data  (rp_rd_data),
        .rp_rd_valid (rp_rd_valid),
        .en_wr_en    (en_wr_en),
        .en_wr_addr  (en_wr_addr),
        .en_wr_data  (en_wr_data),
        .en_rd_en    (en_rd_en),
        .en_rd_addr  (en_rd_addr),
        .en_rd_data  (en_rd_data),
        .en_rd_valid (en_rd_valid)
    );

    // =========================================================================
    // sparsity_scanner
    // =========================================================================
    logic scan_start;
    logic [CSR_ADDR_W-1:0] scan_nnz;

    sparsity_scanner u_scan (
        .clk         (clk),
        .rst_n       (rst_n),
        .start       (scan_start),
        .num_edges   (num_edges),
        .active_n    (active_n),
        .mm_rd_en    (sc_rd_en),
        .mm_rd_addr  (sc_rd_addr),
        .mm_rd_data  (mm_rd_data),
        .mm_rd_valid (mm_rd_valid),
        .rp_wr_en    (rp_wr_en),
        .rp_wr_row   (rp_wr_row),
        .rp_wr_data  (rp_wr_data),
        .en_wr_en    (en_wr_en),
        .en_wr_addr  (en_wr_addr),
        .en_wr_data  (en_wr_data),
        .done        (scan_done),
        .nnz_total   (scan_nnz)
    );

    // =========================================================================
    // sparse_core
    // =========================================================================
    logic [$clog2(N_P)-1:0]   cur_osc;
    logic                      core_start;
    logic signed [XY_W_P-1:0] core_x_new, core_y_new;
    logic                      core_sign_new;
    logic                      core_done;

    sparse_core #(
        .N_P        (N_P),
        .IC_BITS_P  (IC_BITS_P),
        .XY_W_P     (XY_W_P),
        .XY_FRAC_P  (XY_FRAC_P),
        .A_BITS_P   (A_BITS_P),
        .ACCUM_W_P  (ACCUM_W_P),
        .PROD_W_P   (PROD_W_P)
    ) u_core (
        .clk         (clk),
        .rst_n       (rst_n),
        .start       (core_start),
        .osc_idx     (cur_osc),
        .x_i         (x_reg[cur_osc]),
        .y_i         (y_reg[cur_osc]),
        .a_m         (a_m),
        .dt_fp       (dt_fp),
        .c0_fp       (c0_fp),
        .osc_signs   (signs_out),
        .rp_rd_en    (rp_rd_en),
        .rp_rd_row   (rp_rd_row),
        .rp_rd_data  (rp_rd_data),
        .rp_rd_valid (rp_rd_valid),
        .en_rd_en    (en_rd_en),
        .en_rd_addr  (en_rd_addr),
        .en_rd_data  (en_rd_data),
        .en_rd_valid (en_rd_valid),
        .mm_rd_en    (core_rd_en),
        .mm_rd_addr  (core_rd_addr),
        .mm_rd_data  (mm_rd_data),
        .mm_rd_valid (mm_rd_valid),
        .x_i_new     (core_x_new),
        .y_i_new     (core_y_new),
        .sign_xi_new (core_sign_new),
        .done        (core_done)
    );

    // =========================================================================
    // Controller FSM
    // Identical structure to original dsb_array FSM, with two additions:
    //   BOOT state: fires scanner once when run is first asserted
    //   x/y init write path kept exactly as original
    // =========================================================================
    typedef enum logic [2:0] {
        ARR_IDLE  = 3'd0,
        ARR_BOOT  = 3'd1,   // waiting for scan_done
        ARR_LOAD  = 3'd2,   // same as original CTRL_LOAD
        ARR_RUN   = 3'd3,   // same as original CTRL_RUN
        ARR_WAIT  = 3'd4    // same as original CTRL_IDLE (wait for core_done)
    } arr_state_t;

    arr_state_t ctrl;
    logic running;

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            running       <= 1'b0;
            m             <= '0;
            cur_osc       <= '0;
            step_done     <= 1'b0;
            schedule_done <= 1'b0;
            core_start    <= 1'b0;
            scan_start    <= 1'b0;
            ctrl          <= ARR_IDLE;
            for (int i = 0; i < N_P; i++) begin
                // Default small non-zero seed for ALL N_P oscillator slots.
                // Only the first active_n of these participate in the sweep
                // (see ARR_WAIT below); unused slots are simply never touched
                // by the FSM and never referenced by any CSR entry, since the
                // edge list only contains edges among nodes 0..active_n-1.
                x_reg[i] <= XY_W_P'((i + 1) << (XY_FRAC_P - 4));
                y_reg[i] <= '0;
            end
        end else begin
            step_done  <= 1'b0;
            core_start <= 1'b0;
            scan_start <= 1'b0;

            // ── x/y init write (identical to original) ────────────────────────
            if (wr_xy_en && !running)
                begin
                    x_reg[wr_xy_idx] <= wr_x;
                    y_reg[wr_xy_idx] <= wr_y;
                end

            unique case (ctrl)

                // ── Idle: wait for run ────────────────────────────────────────
                ARR_IDLE: begin
                    if (run && !running && !schedule_done) begin
                        running    <= 1'b1;
                        scan_start <= 1'b1;   // kick sparsity scanner
                        ctrl       <= ARR_BOOT;
                    end
                end

                // ── Boot: wait for scanner to finish ──────────────────────────
                ARR_BOOT: begin
                    if (scan_done) begin
                        m       <= '0;
                        cur_osc <= '0;
                        ctrl    <= ARR_LOAD;
                    end
                end

                // ── Load: identical to original CTRL_LOAD ─────────────────────
                // In the sparse version there is no wwl/wbl push — the core
                // reads J from main_memory itself.  We just advance to RUN.
                ARR_LOAD: begin
                    ctrl <= ARR_RUN;
                end

                // ── Run: fire core_start for one cycle ────────────────────────
                ARR_RUN: begin
                    core_start <= 1'b1;
                    ctrl       <= ARR_WAIT;
                end

                // ── Wait: hold until core signals done ────────────────────────
                ARR_WAIT: begin
                    if (core_done) begin
                        x_reg[cur_osc] <= core_x_new;
                        y_reg[cur_osc] <= core_y_new;

                        if (cur_osc == $clog2(N_P)'(active_n - 1)) begin
                            cur_osc   <= '0;
                            step_done <= 1'b1;
                            if (m == Nstep - 1) begin
                                schedule_done <= 1'b1;
                                running       <= 1'b0;
                                ctrl          <= ARR_IDLE;
                            end else begin
                                m    <= m + 1;
                                ctrl <= ARR_LOAD;
                            end
                        end else begin
                            cur_osc <= cur_osc + 1;
                            ctrl    <= ARR_LOAD;
                        end
                    end
                end

                default: ctrl <= ARR_IDLE;

            endcase
        end
    end

endmodule

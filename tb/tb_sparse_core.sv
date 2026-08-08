// tb_sparse_core.sv — corrected with proper fixed‑point reference model
`timescale 1ns/1ps

module tb_sparse_core;

    import dsb_pkg::*;
    import dsb_sparse_pkg::*;

    // ── Local parameters ──────────────────────────────────────────────
    localparam int N_TEST          = 4;
    localparam int K_MAX_P         = 2;
    localparam int IC_BITS_P       = 4;
    localparam int COL_IDX_WP      = $clog2(N_TEST);
    localparam int MAX_EDGES_TEST  = 4;
    localparam int MAX_NNZ_TEST    = 2 * MAX_EDGES_TEST;
    localparam int MAIN_ADDR_WP    = $clog2(MAX_EDGES_TEST + 2*N_TEST);
    localparam int CSR_ADDR_W_P    = $clog2(MAX_NNZ_TEST + 1);
    localparam int CSR_ENTRY_WP    = COL_IDX_WP + MAIN_ADDR_WP;
    localparam int ROW_COUNT_WP    = $clog2(MAX_NNZ_TEST + 1);
    localparam int ROW_PTR_W_P     = CSR_ADDR_W_P + ROW_COUNT_WP;
    localparam int MAIN_WORD_WP    = (COL_IDX_WP + COL_IDX_WP + IC_BITS_P > XY_W) ?
                                      (COL_IDX_WP + COL_IDX_WP + IC_BITS_P) : XY_W;
    localparam int MAIN_MEM_WORDS  = MAX_EDGES_TEST + 2*N_TEST;
    localparam int K_CNT_W_P       = $clog2(K_MAX_P + 1);

    // Clock
    logic clk = 0;
    always #5 clk = ~clk;
    logic rst_n = 1;

    // ── DUT signals ──────────────────────────────────────────────────
    logic                         start;
    logic [COL_IDX_WP-1:0]        osc_idx;
    logic signed [XY_W-1:0]       x_i;
    logic signed [XY_W-1:0]       y_i;
    logic [A_BITS-1:0]            a_m;
    logic [XY_FRAC-1:0]           dt_fp;
    logic [XY_FRAC-1:0]           c0_fp;
    logic [N_TEST-1:0]            osc_signs;
    logic                         rp_rd_en;
    logic [COL_IDX_WP-1:0]        rp_rd_row;
    logic [ROW_PTR_W_P-1:0]       rp_rd_data;
    logic                         rp_rd_valid;
    logic                         en_rd_en;
    logic [CSR_ADDR_W_P-1:0]      en_rd_addr;
    logic [CSR_ENTRY_WP-1:0]      en_rd_data;
    logic                         en_rd_valid;
    logic                         mm_rd_en;
    logic [MAIN_ADDR_WP-1:0]      mm_rd_addr;
    logic [MAIN_WORD_WP-1:0]      mm_rd_data;
    logic                         mm_rd_valid;
    logic signed [XY_W-1:0]       x_i_new;
    logic signed [XY_W-1:0]       y_i_new;
    logic                         sign_xi_new;
    logic                         done;

    // ── DUT ──────────────────────────────────────────────────────────
    sparse_core #(
        .N_P        (N_TEST),
        .IC_BITS_P  (IC_BITS_P),
        .XY_W_P     (XY_W),
        .XY_FRAC_P  (XY_FRAC),
        .A_BITS_P   (A_BITS),
        .ACCUM_W_P  (ACCUM_W),
        .PROD_W_P   (PROD_W),
        .K_MAX_P    (K_MAX_P),
        .COL_IDX_WP (COL_IDX_WP),
        .MAIN_AW    (MAIN_ADDR_WP),
        .MAIN_WW    (MAIN_WORD_WP),
        .CSR_AW     (CSR_ADDR_W_P),
        .CSR_EW     (CSR_ENTRY_WP),
        .ROW_PTR_WP (ROW_PTR_W_P),
        .ROW_CNT_WP (ROW_COUNT_WP)
    ) u_core (
        .clk         (clk),
        .rst_n       (rst_n),
        .start       (start),
        .osc_idx     (osc_idx),
        .x_i         (x_i),
        .y_i         (y_i),
        .a_m         (a_m),
        .dt_fp       (dt_fp),
        .c0_fp       (c0_fp),
        .osc_signs   (osc_signs),
        .rp_rd_en    (rp_rd_en),
        .rp_rd_row   (rp_rd_row),
        .rp_rd_data  (rp_rd_data),
        .rp_rd_valid (rp_rd_valid),
        .en_rd_en    (en_rd_en),
        .en_rd_addr  (en_rd_addr),
        .en_rd_data  (en_rd_data),
        .en_rd_valid (en_rd_valid),
        .mm_rd_en    (mm_rd_en),
        .mm_rd_addr  (mm_rd_addr),
        .mm_rd_data  (mm_rd_data),
        .mm_rd_valid (mm_rd_valid),
        .x_i_new     (x_i_new),
        .y_i_new     (y_i_new),
        .sign_xi_new (sign_xi_new),
        .done        (done)
    );

    // ── External main_memory and csr_index_store ────────────────────
    logic mm_wr_en;
    logic [MAIN_ADDR_WP-1:0] mm_wr_addr;
    logic [MAIN_WORD_WP-1:0] mm_wr_data;

    main_memory #(
        .WORDS_P  (MAIN_MEM_WORDS),
        .WORD_W_P (MAIN_WORD_WP),
        .ADDR_W_P (MAIN_ADDR_WP)
    ) u_mm (
        .clk     (clk),
        .rst_n   (rst_n),
        .wr_en   (mm_wr_en),
        .wr_addr (mm_wr_addr),
        .wr_data (mm_wr_data),
        .rd_en   (mm_rd_en),
        .rd_addr (mm_rd_addr),
        .rd_data (mm_rd_data),
        .rd_valid(mm_rd_valid)
    );

    logic rp_wr_en;
    logic [COL_IDX_WP-1:0] rp_wr_row;
    logic [ROW_PTR_W_P-1:0] rp_wr_data;
    logic en_wr_en;
    logic [CSR_ADDR_W_P-1:0] en_wr_addr;
    logic [CSR_ENTRY_WP-1:0] en_wr_data;

    csr_index_store #(
        .N_P          (N_TEST),
        .MAX_NNZ_P    (MAX_NNZ_TEST),
        .COL_IDX_W_P  (COL_IDX_WP),
        .MAIN_ADDR_WP (MAIN_ADDR_WP),
        .CSR_ADDR_W_P (CSR_ADDR_W_P),
        .CSR_ENTRY_WP (CSR_ENTRY_WP),
        .ROW_PTR_W_P  (ROW_PTR_W_P),
        .ROW_COUNT_WP (ROW_COUNT_WP)
    ) u_csr (
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

    // ── Pre‑load tasks ──────────────────────────────────────────────
    task mm_write(input int addr, input logic [MAIN_WORD_WP-1:0] data);
        @(posedge clk);
        mm_wr_en = 1; mm_wr_addr = addr; mm_wr_data = data;
        @(posedge clk);
        mm_wr_en = 0;
    endtask

    task rp_write(input int row, input int start, input int count);
        logic [ROW_PTR_W_P-1:0] data = {CSR_ADDR_W_P'(start), ROW_COUNT_WP'(count)};
        @(posedge clk);
        rp_wr_en = 1; rp_wr_row = row; rp_wr_data = data;
        @(posedge clk);
        rp_wr_en = 0;
    endtask

    task en_write(input int addr, input int col, input int maddr);
        logic [CSR_ENTRY_WP-1:0] data = {COL_IDX_WP'(col), MAIN_ADDR_WP'(maddr)};
        @(posedge clk);
        en_wr_en = 1; en_wr_addr = addr; en_wr_data = data;
        @(posedge clk);
        en_wr_en = 0;
    endtask

    function automatic logic [MAIN_WORD_WP-1:0] pack_edge(input int r, c, w);
        return {COL_IDX_WP'(r), COL_IDX_WP'(c), IC_BITS_P'(signed'(w))};
    endfunction

    // ── Reference model for Jx_i ────────────────────────────────────
    function automatic logic signed [ACCUM_W-1:0] compute_Jx(
        input logic [N_TEST-1:0] signs
    );
        logic signed [ACCUM_W-1:0] Jx = 0;
        // Hardcoded edges: (0,1)=5, (0,2)=-3, (0,3)=2
        Jx += 5 * (signs[1] ? 1 : -1);
        Jx += (-3) * (signs[2] ? 1 : -1);
        Jx += 2 * (signs[3] ? 1 : -1);
        return Jx;
    endfunction

    // ── Correct reference update_unit (matches RTL fixed‑point) ─────
    function automatic void ref_update(
        input logic signed [XY_W-1:0] x,
        input logic signed [XY_W-1:0] y,
        input logic signed [ACCUM_W-1:0] Jx,
        input logic [A_BITS-1:0] a,
        input logic [XY_FRAC-1:0] dt,
        input logic [XY_FRAC-1:0] c0,
        output logic signed [XY_W-1:0] x_new,
        output logic signed [XY_W-1:0] y_new,
        output logic sign_new
    );
        // Stage A: a*x
        logic signed [PROD_W-1:0] ax_full;
        logic signed [XY_W-1:0] ax_shifted;
        ax_full = $signed({{(PROD_W-XY_W){x[XY_W-1]}}, x}) * $signed({1'b0, a});
        ax_shifted = XY_W'(ax_full >>> (A_BITS-1));

        // Stage B: c0*Jx
        logic signed [ACCUM_W+XY_FRAC-1:0] c0jx_full;
        c0jx_full = Jx * $signed({1'b0, c0});

        // Stage C: Δy
        logic signed [XY_W-1:0] neg_ax = -ax_shifted;
        logic signed [ACCUM_W+XY_FRAC-1:0] force_wide;
        force_wide = $signed({{(ACCUM_W+XY_FRAC-XY_W){neg_ax[XY_W-1]}}, neg_ax})
                     + $signed(c0jx_full);
        logic signed [ACCUM_W+2*XY_FRAC-1:0] dy_full;
        dy_full = force_wide * $signed({1'b0, dt});
        logic signed [XY_W-1:0] delta_y;
        delta_y = XY_W'(dy_full >>> XY_FRAC);

        // Stage D: y_pre
        logic signed [XY_W-1:0] y_pre = y + delta_y;

        // Stage E: x_pre
        logic signed [2*XY_W-1:0] dx_full;
        dx_full = y_pre * $signed({1'b0, dt});
        logic signed [XY_W-1:0] x_pre = x + XY_W'(dx_full >>> XY_FRAC);

        // Stage F: wall
        logic wall_hit_pos = (x_pre >  $signed(ONE_FP));
        logic wall_hit_neg = (x_pre < -$signed(ONE_FP));
        logic wall_hit = wall_hit_pos | wall_hit_neg;
        if (wall_hit) begin
            x_new = wall_hit_pos ?  ONE_FP : -ONE_FP;
            y_new = '0;
        end else begin
            x_new = x_pre;
            y_new = y_pre;
        end
        sign_new = ~x_new[XY_W-1];
    endfunction

    // ── Test control ──────────────────────────────────────────────────
    int error_count = 0;

    task check_equal(input string msg, input logic signed [XY_W-1:0] a,
                     input logic signed [XY_W-1:0] b);
        if (a !== b) begin
            $error("%s: expected %0d, got %0d", msg, b, a);
            error_count++;
        end else begin
            $display("  %s: OK", msg);
        end
    endtask

    task check_equal_sign(input string msg, input logic a, input logic b);
        if (a !== b) begin
            $error("%s: expected %b, got %b", msg, b, a);
            error_count++;
        end else begin
            $display("  %s: OK", msg);
        end
    endtask

    // ── Test sequence ──────────────────────────────────────────────────
    initial begin
        // Initialize
        mm_wr_en   = 0;
        rp_wr_en   = 0;
        en_wr_en   = 0;
        start = 0;
        osc_idx = 0;
        x_i = 50 * ONE_FP;    // 50.0 in Q1.14
        y_i = 10 * ONE_FP;    // 10.0
        a_m = 8192;           // 0.5
        dt_fp = 164;          // 0.01
        c0_fp = 8192;         // 0.5
        osc_signs = 4'b0010;  // node0=+1, node1=+1, node2=-1, node3=+1

        repeat (4) @(posedge clk);
        rst_n = 1;
        repeat (2) @(posedge clk);
        $display("=== Starting sparse_core test ===");

        // ── Pre‑load main memory ──────────────────────────────────────────
        $display("Pre‑load main memory");
        mm_write(0, pack_edge(0,1,5));
        mm_write(1, pack_edge(0,2,-3));
        mm_write(2, pack_edge(0,3,2));

        // ── Pre‑load CSR store ──────────────────────────────────────────
        $display("Pre‑load CSR store");
        rp_write(0, 0, 3);
        en_write(0, 1, 0);
        en_write(1, 2, 1);
        en_write(2, 3, 2);

        // ── Compute expected result ──────────────────────────────────────
        logic signed [ACCUM_W-1:0] Jx_correct;
        Jx_correct = compute_Jx(osc_signs);  // 6
        logic signed [XY_W-1:0] x_exp, y_exp;
        logic sign_exp;
        ref_update(x_i, y_i, Jx_correct, a_m, dt_fp, c0_fp, x_exp, y_exp, sign_exp);
        $display("Expected Jx_i = %0d, x_new=%0d, y_new=%0d, sign=%b",
                 Jx_correct, x_exp, y_exp, sign_exp);

        // ── Run the core for oscillator 0 ──────────────────────────────
        $display("Starting core...");
        @(posedge clk);
        start = 1;
        @(posedge clk);
        start = 0;

        // Wait for done (timeout)
        repeat (500) @(posedge clk) if (done) break;
        if (!done) $fatal("Timeout waiting for done");

        // ── Check outputs ──────────────────────────────────────────────
        $display("Checking outputs...");
        check_equal("x_i_new", x_i_new, x_exp);
        check_equal("y_i_new", y_i_new, y_exp);
        check_equal_sign("sign_xi_new", sign_xi_new, sign_exp);

        // ── Summary ─────────────────────────────────────────────────────
        if (error_count == 0) begin
            $display("=== All tests PASSED ===");
        end else begin
            $error("=== %0d tests FAILED ===", error_count);
            $finish(1);
        end
        $finish;
    end

    // VCD dump
    initial begin
        $dumpfile("tb_sparse_core.vcd");
        $dumpvars(0, tb_sparse_core);
    end

endmodule
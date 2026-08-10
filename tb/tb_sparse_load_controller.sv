// tb_sparse_load_controller.sv
`timescale 1ns/1ps

module tb_sparse_load_controller;

    import dsb_pkg::*;
    import dsb_sparse_pkg::*;

    // ── Local parameters (small test system) ──────────────────────────
    localparam int N_TEST          = 4;
    localparam int K_MAX_P         = 2;            // force chunking
    localparam int IC_BITS_P       = 4;
    localparam int COL_IDX_W_P     = $clog2(N_TEST);
    localparam int MAX_EDGES_TEST  = 4;
    localparam int MAX_NNZ_TEST    = 2 * MAX_EDGES_TEST;
    localparam int MAIN_ADDR_WP    = $clog2(MAX_EDGES_TEST + 2*N_TEST);
    localparam int CSR_ADDR_W_P    = $clog2(MAX_NNZ_TEST + 1);
    localparam int CSR_ENTRY_WP    = COL_IDX_W_P + MAIN_ADDR_WP;
    localparam int ROW_COUNT_WP    = $clog2(MAX_NNZ_TEST + 1);
    localparam int ROW_PTR_W_P     = CSR_ADDR_W_P + ROW_COUNT_WP;
    localparam int MAIN_WORD_WP    = (COL_IDX_W_P+COL_IDX_W_P+IC_BITS_P > XY_W) ?
                                      (COL_IDX_W_P+COL_IDX_W_P+IC_BITS_P) : XY_W;
    localparam int MAIN_MEM_WORDS  = MAX_EDGES_TEST + 2*N_TEST;
    localparam int K_CNT_W_P       = $clog2(K_MAX_P+1);
    localparam int ACCUM_W_P       = 24;

    logic clk = 0;
    always #5 clk = ~clk;
    logic rst_n = 1;

    // ── DUT signals ──────────────────────────────────────────────────
    logic                        slc_start;
    logic [COL_IDX_W_P-1:0]      osc_idx;
    logic                        slc_done;
    logic [N_TEST-1:0]           osc_signs;

    // CSR read ports
    logic                        rp_rd_en;
    logic [COL_IDX_W_P-1:0]      rp_rd_row;
    logic [ROW_PTR_W_P-1:0]      rp_rd_data;
    logic                        rp_rd_valid;
    logic                        en_rd_en;
    logic [CSR_ADDR_W_P-1:0]     en_rd_addr;
    logic [CSR_ENTRY_WP-1:0]     en_rd_data;
    logic                        en_rd_valid;

    // Main memory read ports
    logic                        mm_rd_en;
    logic [MAIN_ADDR_WP-1:0]     mm_rd_addr;
    logic [MAIN_WORD_WP-1:0]     mm_rd_data;
    logic                        mm_rd_valid;

    // SCM write
    logic                        scm_wr_en;
    logic [$clog2(K_MAX_P)-1:0]  scm_wr_slot;
    logic [COL_IDX_W_P-1:0]      scm_wr_col_idx;
    logic [IC_BITS_P-1:0]        scm_wr_J;
    logic                        scm_wr_sign_j;
    logic                        scm_precharge;
    logic [K_MAX_P-1:0]          scm_rwl;

    // Dotprod control (SLC outputs)
    logic                        ph2_start;
    logic                        ph2_accumulate;
    logic                        ph2_last_chunk;
    logic [K_CNT_W_P-1:0]        ph2_valid_count;
    logic                        ph2_sign_xi;
    logic                        ph2_done;

    // ── Instantiate real memory and compute modules ──────────────────
    // Main memory (with write ports for preload)
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

    // CSR store (with write ports for preload)
    logic rp_wr_en;
    logic [COL_IDX_W_P-1:0] rp_wr_row;
    logic [ROW_PTR_W_P-1:0] rp_wr_data;
    logic en_wr_en;
    logic [CSR_ADDR_W_P-1:0] en_wr_addr;
    logic [CSR_ENTRY_WP-1:0] en_wr_data;

    csr_index_store #(
        .N_P          (N_TEST),
        .MAX_NNZ_P    (MAX_NNZ_TEST),
        .COL_IDX_W_P  (COL_IDX_W_P),
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

    // Sparse compute memory
    logic [K_MAX_P*IC_BITS_P-1:0] rbl_J;
    logic [K_MAX_P-1:0]           rbl_signs;

    sparse_compute_memory #(
        .K_MAX_P  (K_MAX_P),
        .IC_BITS_P(IC_BITS_P),
        .COL_W_P  (COL_IDX_W_P)
    ) u_scm (
        .clk         (clk),
        .wr_en       (scm_wr_en),
        .wr_slot     (scm_wr_slot),
        .wr_col_idx  (scm_wr_col_idx),
        .wr_J        (scm_wr_J),
        .wr_sign_j   (scm_wr_sign_j),
        .rd_col_en   (1'b0),
        .rd_slot     ('0),
        .rd_col_idx  (),
        .rd_col_valid(),
        .precharge   (scm_precharge),
        .rwl         (scm_rwl),
        .rbl_J       (rbl_J),
        .rbl_signs   (rbl_signs)
    );

    // Dotprod
    logic signed [ACCUM_W_P-1:0] Jx_i_out;

    dotprod_phase2_sparse #(
        .K_MAX_P   (K_MAX_P),
        .IC_BITS_P (IC_BITS_P),
        .ACCUM_W_P (ACCUM_W_P)
    ) u_ph2 (
        .clk         (clk),
        .start       (ph2_start),
        .done        (ph2_done),
        .accumulate  (ph2_accumulate),
        .last_chunk  (ph2_last_chunk),
        .valid_count (ph2_valid_count),
        .sign_xi     (ph2_sign_xi),
        .xnor_J      (rbl_J),
        .sign_eq     (rbl_signs),
        .Jx_i        (Jx_i_out)
    );

    // SLC DUT
    sparse_load_controller #(
        .N_P          (N_TEST),
        .K_MAX_P      (K_MAX_P),
        .IC_BITS_P    (IC_BITS_P),
        .COL_IDX_W_P  (COL_IDX_W_P),
        .MAIN_ADDR_WP (MAIN_ADDR_WP),
        .MAIN_WORD_WP (MAIN_WORD_WP),
        .CSR_ADDR_W_P (CSR_ADDR_W_P),
        .CSR_ENTRY_WP (CSR_ENTRY_WP),
        .ROW_PTR_W_P  (ROW_PTR_W_P),
        .ROW_COUNT_WP (ROW_COUNT_WP),
        .K_CNT_W_P    (K_CNT_W_P),
        .XY_W_P       (XY_W)
    ) u_slc (
        .clk            (clk),
        .rst_n          (rst_n),
        .start          (slc_start),
        .osc_idx        (osc_idx),
        .done           (slc_done),
        .osc_signs      (osc_signs),
        .rp_rd_en       (rp_rd_en),
        .rp_rd_row      (rp_rd_row),
        .rp_rd_data     (rp_rd_data),
        .rp_rd_valid    (rp_rd_valid),
        .en_rd_en       (en_rd_en),
        .en_rd_addr     (en_rd_addr),
        .en_rd_data     (en_rd_data),
        .en_rd_valid    (en_rd_valid),
        .mm_rd_en       (mm_rd_en),
        .mm_rd_addr     (mm_rd_addr),
        .mm_rd_data     (mm_rd_data),
        .mm_rd_valid    (mm_rd_valid),
        .scm_wr_en      (scm_wr_en),
        .scm_wr_slot    (scm_wr_slot),
        .scm_wr_col_idx (scm_wr_col_idx),
        .scm_wr_J       (scm_wr_J),
        .scm_wr_sign_j  (scm_wr_sign_j),
        .scm_precharge  (scm_precharge),
        .scm_rwl        (scm_rwl),
        .ph2_start      (ph2_start),
        .ph2_accumulate (ph2_accumulate),
        .ph2_last_chunk (ph2_last_chunk),
        .ph2_valid_count(ph2_valid_count),
        .ph2_sign_xi    (ph2_sign_xi),
        .ph2_done       (ph2_done)
    );

    // ── Helper tasks ──────────────────────────────────────────────────
    task mm_write(input int addr, input logic [MAIN_WORD_WP-1:0] data);
        @(posedge clk);
        mm_wr_en = 1; mm_wr_addr = addr; mm_wr_data = data;
        @(posedge clk);
        mm_wr_en = 0;
    endtask

    function automatic logic [MAIN_WORD_WP-1:0] pack_edge(input int r, c, w);
        return {COL_IDX_W_P'(r), COL_IDX_W_P'(c), IC_BITS_P'(signed'(w))};
    endfunction

    task rp_write(input int row, input int start, input int count);
        automatic logic [ROW_PTR_W_P-1:0] data = {CSR_ADDR_W_P'(start), ROW_COUNT_WP'(count)};
        @(posedge clk);
        rp_wr_en = 1; rp_wr_row = row; rp_wr_data = data;
        @(posedge clk);
        rp_wr_en = 0;
    endtask

    task en_write(input int addr, input int col, input int maddr);
        automatic logic [CSR_ENTRY_WP-1:0] data = {COL_IDX_W_P'(col), MAIN_ADDR_WP'(maddr)};
        @(posedge clk);
        en_wr_en = 1; en_wr_addr = addr; en_wr_data = data;
        @(posedge clk);
        en_wr_en = 0;
    endtask

    // ── Test control ──────────────────────────────────────────────────
    int error_count = 0;

    task check_equal_val(input string msg, input logic a, input logic b);
        if (a !== b) begin
            $error("%s: expected %b, got %b", msg, b, a);
            error_count++;
        end else begin
            $display("  %s: OK", msg);
        end
    endtask

    task check_equal_int(input string msg, input int a, input int b);
        if (a !== b) begin
            $error("%s: expected %0d, got %0d", msg, b, a);
            error_count++;
        end else begin
            $display("  %s: OK", msg);
        end
    endtask

    // ── Test sequence ──────────────────────────────────────────────────
    initial begin
        // Initialize
        mm_wr_en = 0; rp_wr_en = 0; en_wr_en = 0;
        osc_idx = 0;
        osc_signs = 4'b0010; // node0=+1, node1=+1, node2=0, node3=+1
        slc_start = 0;

        repeat (4) @(posedge clk);
        rst_n = 1;
        repeat (2) @(posedge clk);
        $display("=== sparse_load_controller test (RWL=1, valid_count=chunk_valid) ===");

        // ── Preload main memory ──────────────────────────────────────
        // Row 0 has 3 edges: (0,1) J=5, (0,2) J=-3, (0,3) J=2
        mm_write(0, pack_edge(0,1,5));
        mm_write(1, pack_edge(0,2,-3));
        mm_write(2, pack_edge(0,3,2));

        // ── Preload CSR store ─────────────────────────────────────────
        rp_write(0, 0, 3);
        en_write(0, 1, 0);
        en_write(1, 2, 1);
        en_write(2, 3, 2);

        // ── Start SLC for row 0 (osc_idx=0) ────────────────────────────
        @(posedge clk);
        slc_start = 1;
        @(posedge clk);
        slc_start = 0;

        // Wait for done (timeout)
        repeat (100) @(posedge clk) if (slc_done) break;
        if (!slc_done) $fatal("Timeout waiting for slc_done");

        // ── Verify ph2 control signals ─────────────────────────────────
        // We will count the number of ph2_start pulses and check accumulate, valid_count.
        // Since K_MAX=2, row0 has 3 entries -> two chunks: chunk0 (2 entries), chunk1 (1 entry)
        // We'll use a monitor (below) to record the parameters and check them.

        // ── After run, check final Jx_i ──────────────────────────────
        // Expected with RWL=1 and case table:
        // For each term, with RWL=1, the SCM outputs XNOR(J, 1) = J and sign_eq = XNOR(sign_j, 1) = sign_j.
        // The dotprod case table will then compute J * sign_j correctly.
        // Row0: J01=5, sign_j1=+1 -> product 5; J02=-3, sign_j2=-1 -> product 3; J03=2, sign_j3=+1 -> product 2; sum=10.
        // However the SLC only loads all 3 entries and triggers dotprod twice (chunk0 and chunk1).
        // The final Jx_i should be 10.
        check_equal_int("Final Jx_i", $signed(Jx_i_out), 10);

        // ── Summary ──────────────────────────────────────────────────
        if (error_count == 0) begin
            $display("=== All tests PASSED ===");
        end else begin
            $error("=== %0d tests FAILED ===", error_count);
            $finish(1);
        end
        $finish;
    end

    // ── Monitor for ph2_start parameters ────────────────────────────
    integer chunk_count = 0;
    logic [K_CNT_W_P-1:0] expected_valid[2] = '{2, 1}; // chunk0=2, chunk1=1
    logic expected_accum[2] = '{0, 1};
    logic expected_last[2] = '{0, 1};

    always @(posedge clk) begin
        if (ph2_start) begin
            $display("[%0t] ph2_start #%0d: accum=%b, last=%b, valid=%0d, sign_xi=%b",
                     $time, chunk_count+1, ph2_accumulate, ph2_last_chunk,
                     ph2_valid_count, ph2_sign_xi);
            // Check scm_rwl = all 1's
            if (scm_rwl !== {K_MAX_P{1'b1}}) begin
                $error("scm_rwl not all 1's, got %b", scm_rwl);
                error_count++;
            end else begin
                $display("  scm_rwl = all 1's OK");
            end
            // Check valid_count and accumulate
            if (chunk_count < 2) begin
                if (ph2_valid_count !== expected_valid[chunk_count]) begin
                    $error("chunk%d valid_count: expected %0d, got %0d",
                           chunk_count, expected_valid[chunk_count], ph2_valid_count);
                    error_count++;
                end
                if (ph2_accumulate !== expected_accum[chunk_count]) begin
                    $error("chunk%d accumulate: expected %b, got %b",
                           chunk_count, expected_accum[chunk_count], ph2_accumulate);
                    error_count++;
                end
                if (ph2_last_chunk !== expected_last[chunk_count]) begin
                    $error("chunk%d last_chunk: expected %b, got %b",
                           chunk_count, expected_last[chunk_count], ph2_last_chunk);
                    error_count++;
                end
            end
            chunk_count++;
        end
    end

    initial begin
        $dumpfile("tb_sparse_load_controller.vcd");
        $dumpvars(0, tb_sparse_load_controller);
    end

endmodule
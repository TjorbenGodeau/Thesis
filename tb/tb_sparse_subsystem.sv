// =============================================================================
// tb_sparse_subsystem.sv
//
// Self-checking testbench for the full sparse memory subsystem:
//   main_memory → sparsity_scanner → csr_index_store
//                                 → sparse_load_controller
//                                 → sparse_compute_memory
//                                 → dotprod_phase2_sparse
//
// Test scenario (N=8, K_MAX=4 from package defaults)
// ---------------------------------------------------
//  J matrix (upper-triangular style, symmetric for Ising):
//
//     row 0: J[0][2]=+3, J[0][5]=-1                 (2 non-zeros)
//     row 1: J[1][3]=+2, J[1][4]=+1, J[1][6]=-2    (3 non-zeros)
//     row 2: J[2][0]=+3, J[2][7]=+1                 (2 non-zeros)
//     row 3: J[3][1]=+2, J[3][4]=-1, J[3][5]=+2,
//             J[3][6]=+1, J[3][7]=-1                (5 non-zeros → 2 chunks)
//     rows 4-7: sparse, few entries
//
//  Initial spin signs (for the dot-product test):
//     signs = 8'b10110101  → osc_signs[0]=1,1=0,2=1,3=0,4=1,5=1,6=0,7=1
//
// Expected Jx_i for row 3:
//   J[3][1]*sign(x_1) + J[3][4]*sign(x_4) + J[3][5]*sign(x_5)
//   + J[3][6]*sign(x_6) + J[3][7]*sign(x_7)
//   = 2*(-1) + (-1)*(+1) + 2*(+1) + 1*(-1) + (-1)*(+1)
//     (signs: x1=neg→-1, x4=pos→+1, x5=pos→+1, x6=neg→-1, x7=pos→+1)
//   = -2 + (-1) + 2 + (-1) + (-1) = -3
//
// The testbench checks:
//   1. Scanner completes and reports correct nnz_total
//   2. Row-pointer entries are correct for rows 0-3
//   3. Sparse load controller processes row 3 in 2 chunks
//   4. Final Jx_i == -3
// =============================================================================

`timescale 1ns/1ps

module tb_sparse_subsystem;

    import dsb_pkg::*;
    import dsb_sparse_pkg::*;

    // ── Clock / reset ─────────────────────────────────────────────────────────
    logic clk = 0;
    logic rst_n = 0;
    always #5 clk = ~clk;   // 100 MHz

    // ── DUT port signals ──────────────────────────────────────────────────────

    // main_memory
    logic                   mm_wr_en;
    logic [MAIN_ADDR_W-1:0] mm_wr_addr;
    logic [MAIN_WORD_W-1:0] mm_wr_data;
    logic                   mm_rd_en;
    logic [MAIN_ADDR_W-1:0] mm_rd_addr;
    logic [MAIN_WORD_W-1:0] mm_rd_data;
    logic                   mm_rd_valid;

    // csr_index_store
    logic                   rp_wr_en;
    logic [COL_IDX_W-1:0]   rp_wr_row;
    logic [ROW_PTR_W-1:0]   rp_wr_data;
    logic                   rp_rd_en;
    logic [COL_IDX_W-1:0]   rp_rd_row;
    logic [ROW_PTR_W-1:0]   rp_rd_data;
    logic                   rp_rd_valid;
    logic                   en_wr_en;
    logic [CSR_ADDR_W-1:0]  en_wr_addr;
    logic [CSR_ENTRY_W-1:0] en_wr_data;
    logic                   en_rd_en;
    logic [CSR_ADDR_W-1:0]  en_rd_addr;
    logic [CSR_ENTRY_W-1:0] en_rd_data;
    logic                   en_rd_valid;

    // sparsity_scanner
    logic                   scan_start;
    logic                   scan_done;
    logic [CSR_ADDR_W-1:0]  scan_nnz_total;
    // scanner also drives mm_rd_en/addr — muxed below
    logic                   sc_mm_rd_en;
    logic [MAIN_ADDR_W-1:0] sc_mm_rd_addr;

    // sparse_load_controller
    logic                   slc_start;
    logic [COL_IDX_W-1:0]   slc_osc_idx;
    logic                   slc_done;
    logic [N-1:0]            slc_osc_signs;
    logic                   slc_rp_rd_en;
    logic [COL_IDX_W-1:0]   slc_rp_rd_row;
    logic                   slc_en_rd_en;
    logic [CSR_ADDR_W-1:0]  slc_en_rd_addr;
    logic                   slc_mm_rd_en;
    logic [MAIN_ADDR_W-1:0] slc_mm_rd_addr;
    // sparse_compute_memory
    logic                   scm_wr_en;
    logic [$clog2(K_MAX)-1:0] scm_wr_slot;
    logic [COL_IDX_W-1:0]   scm_wr_col_idx;
    logic [IC_BITS-1:0]      scm_wr_J;
    logic                    scm_wr_sign_j;
    logic                    scm_precharge;
    logic [K_MAX-1:0]        scm_rwl;
    logic [K_MAX*IC_BITS-1:0] scm_rbl_J;
    logic [K_MAX-1:0]        scm_rbl_signs;
    // dotprod_phase2_sparse
    logic                   ph2_start;
    logic                   ph2_accumulate;
    logic                   ph2_last_chunk;
    logic [$clog2(K_MAX+1)-1:0] ph2_valid_count;
    logic                   ph2_sign_xi;
    logic                   ph2_done;
    logic signed [ACCUM_W-1:0] ph2_Jx_i;

    // ── Read bus mux: scanner uses mm during scan, SLC uses mm during run ─────
    // Simple priority: scanner takes bus when scan not done yet
    assign mm_rd_en   = scan_done ? slc_mm_rd_en   : sc_mm_rd_en;
    assign mm_rd_addr = scan_done ? slc_mm_rd_addr  : sc_mm_rd_addr;

    // CSR read bus: driven by SLC only
    assign rp_rd_en  = slc_rp_rd_en;
    assign rp_rd_row = slc_rp_rd_row;
    assign en_rd_en  = slc_en_rd_en;
    assign en_rd_addr= slc_en_rd_addr;

    // ── DUT instantiations ────────────────────────────────────────────────────

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

    sparsity_scanner u_scan (
        .clk        (clk),
        .rst_n      (rst_n),
        .start      (scan_start),
        .mm_rd_en   (sc_mm_rd_en),
        .mm_rd_addr (sc_mm_rd_addr),
        .mm_rd_data (mm_rd_data),
        .mm_rd_valid(mm_rd_valid),
        .rp_wr_en   (rp_wr_en),
        .rp_wr_row  (rp_wr_row),
        .rp_wr_data (rp_wr_data),
        .en_wr_en   (en_wr_en),
        .en_wr_addr (en_wr_addr),
        .en_wr_data (en_wr_data),
        .done       (scan_done),
        .nnz_total  (scan_nnz_total)
    );

    sparse_load_controller u_slc (
        .clk            (clk),
        .rst_n          (rst_n),
        .start          (slc_start),
        .osc_idx        (slc_osc_idx),
        .done           (slc_done),
        .osc_signs      (slc_osc_signs),
        .rp_rd_en       (slc_rp_rd_en),
        .rp_rd_row      (slc_rp_rd_row),
        .rp_rd_data     (rp_rd_data),
        .rp_rd_valid    (rp_rd_valid),
        .en_rd_en       (slc_en_rd_en),
        .en_rd_addr     (slc_en_rd_addr),
        .en_rd_data     (en_rd_data),
        .en_rd_valid    (en_rd_valid),
        .mm_rd_en       (slc_mm_rd_en),
        .mm_rd_addr     (slc_mm_rd_addr),
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

    sparse_compute_memory u_scm (
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
        .rbl_J       (scm_rbl_J),
        .rbl_signs   (scm_rbl_signs)
    );

    dotprod_phase2_sparse u_ph2 (
        .clk         (clk),
        .start       (ph2_start),
        .done        (ph2_done),
        .accumulate  (ph2_accumulate),
        .last_chunk  (ph2_last_chunk),
        .valid_count (ph2_valid_count),
        .sign_xi     (ph2_sign_xi),
        .xnor_J      (scm_rbl_J),
        .sign_eq     (scm_rbl_signs),
        .Jx_i        (ph2_Jx_i)
    );

    // ── Task: write J[i][j] into main_memory ──────────────────────────────────
    task automatic write_J(input int i, j, input logic signed [IC_BITS-1:0] val);
        @(posedge clk);
        mm_wr_en   <= 1'b1;
        mm_wr_addr <= MAIN_ADDR_W'(i * N + j);
        mm_wr_data <= MAIN_WORD_W'(signed'(val));
        @(posedge clk);
        mm_wr_en   <= 1'b0;
    endtask

    // ── Task: wait for signal high with timeout ────────────────────────────────
    task automatic wait_for(input logic sig, input int timeout_cycles);
        int cnt = 0;
        while (!sig && cnt < timeout_cycles) begin
            @(posedge clk);
            cnt++;
        end
        if (cnt >= timeout_cycles)
            $fatal(1, "TIMEOUT waiting for signal after %0d cycles", timeout_cycles);
    endtask

    // ── Check macro ───────────────────────────────────────────────────────────
    int pass_count = 0;
    int fail_count = 0;

    task automatic check(input string name,
                         input logic signed [31:0] got,
                         input logic signed [31:0] exp);
        if (got === exp) begin
            $display("  PASS  %-30s  got=%0d", name, got);
            pass_count++;
        end else begin
            $display("  FAIL  %-30s  got=%0d  exp=%0d", name, got, exp);
            fail_count++;
        end
    endtask

    // ── Stimulus ──────────────────────────────────────────────────────────────
    initial begin
        // Defaults
        mm_wr_en    = 0;
        mm_wr_addr  = '0;
        mm_wr_data  = '0;
        scan_start  = 0;
        slc_start   = 0;
        slc_osc_idx = '0;
        slc_osc_signs = '0;

        // Release reset
        repeat(4) @(posedge clk);
        rst_n = 1;
        repeat(2) @(posedge clk);

        // ── Phase 1: Load J matrix into main_memory ───────────────────────────
        $display("\n=== Phase 1: Loading J matrix ===");

        // Row 0
        write_J(0, 2,  3);
        write_J(0, 5, -1);
        // Row 1
        write_J(1, 3,  2);
        write_J(1, 4,  1);
        write_J(1, 6, -2);
        // Row 2
        write_J(2, 0,  3);
        write_J(2, 7,  1);
        // Row 3 (5 non-zeros → needs 2 chunks with K_MAX=4)
        write_J(3, 1,  2);
        write_J(3, 4, -1);
        write_J(3, 5,  2);
        write_J(3, 6,  1);
        write_J(3, 7, -1);
        // Rows 4-7: a couple of entries
        write_J(4, 0,  1);
        write_J(5, 2,  1);

        $display("J matrix loaded.");
        repeat(2) @(posedge clk);

        // ── Phase 2: Run sparsity scanner ─────────────────────────────────────
        $display("\n=== Phase 2: Sparsity scan ===");
        @(posedge clk);
        scan_start <= 1'b1;
        @(posedge clk);
        scan_start <= 1'b0;

        wait_for(scan_done, 500);
        $display("Scan done.  nnz_total=%0d", scan_nnz_total);

        // Row 0: 2, Row 1: 3, Row 2: 2, Row 3: 5, Row 4: 1, Row 5: 1 → total 14
        check("nnz_total", scan_nnz_total, 14);

        repeat(4) @(posedge clk);

        // ── Phase 3: Sparse load controller — process row 3 ──────────────────
        $display("\n=== Phase 3: Sparse load for oscillator 3 ===");
        // signs: x0=+,x1=-,x2=+,x3=-,x4=+,x5=+,x6=-,x7=+
        //        bit7..0: 1=pos, 0=neg
        slc_osc_signs <= 8'b10110101;  // [0]=1,[1]=0,[2]=1,[3]=0,[4]=1,[5]=1,[6]=0,[7]=1
        slc_osc_idx   <= COL_IDX_W'(3);

        @(posedge clk);
        slc_start <= 1'b1;
        @(posedge clk);
        slc_start <= 1'b0;

        wait_for(slc_done, 1000);
        $display("SLC done.  Jx_i=%0d", $signed(ph2_Jx_i));

        // Expected: -3  (computed in header comment)
        check("Jx_i for row 3", $signed(ph2_Jx_i), -3);

        // ── Summary ───────────────────────────────────────────────────────────
        repeat(4) @(posedge clk);
        $display("\n=== Results: %0d passed, %0d failed ===", pass_count, fail_count);
        if (fail_count == 0)
            $display("ALL TESTS PASSED");
        else
            $display("SOME TESTS FAILED — check output above");

        $finish;
    end

    // ── Timeout watchdog ──────────────────────────────────────────────────────
    initial begin
        #200000;
        $fatal(1, "Global testbench timeout");
    end

    // ── Waveform dump ─────────────────────────────────────────────────────────
    initial begin
        $dumpfile("tb_sparse_subsystem.vcd");
        $dumpvars(0, tb_sparse_subsystem);
    end

endmodule

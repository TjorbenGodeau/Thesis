// =============================================================================
// tb_gset_benchmark.sv
//
// Configurable benchmark testbench for running ANY Gset instance through the
// sparse dSB Ising machine (sparse_dsb_array) and comparing the resulting
// Max-Cut value against the literature best-known value.
//
// ─────────────────────────────────────────────────────────────────────────────
// HOW TO RUN A DIFFERENT GSET FILE
// ─────────────────────────────────────────────────────────────────────────────
// No recompilation needed to switch instances -- pass the file path and
// instance id as simulator plusargs:
//
//   vvp sparse_dsb_sim +gset_file=G1.txt +gset_id=1 +nstep=50
//
// Defaults (if no plusargs given): G1.txt, id=1, nstep=30.
//
// ─────────────────────────────────────────────────────────────────────────────
// GSET FILE FORMAT  (as published at web.stanford.edu/~yyye/yyye/Gset/)
// ─────────────────────────────────────────────────────────────────────────────
//   Line 1:        <num_nodes> <num_edges>
//   Lines 2..end:  <node_i> <node_j> <weight>      (1-indexed nodes)
//
// This testbench converts node indices to 0-indexed internally.
//
// ─────────────────────────────────────────────────────────────────────────────
// WHAT THIS TESTBENCH DOES
// ─────────────────────────────────────────────────────────────────────────────
//  1. Parse the Gset file with $fopen / $fscanf directly in SystemVerilog
//     (no external preprocessing script needed).
//  2. Load every edge into main_memory via sparse_dsb_array's mm_wr_* port,
//     packed as {row, col, weight} edge words.
//  3. Initialise x with small random +-perturbations (standard dSB practice)
//     via the wr_xy_en port.
//  4. Set active_n = num_nodes, num_edges = num_edges parsed from the file.
//  5. Assert run=1; the array auto-triggers the edge-list scanner, then
//     sweeps Nstep dSB steps.
//  6. On schedule_done, read signs_out, compute the Max-Cut value directly
//     from the edge list (cut = sum of edge weights where the two endpoints
//     have DIFFERENT signs), and compare to gset_bestknown_pkg.
//  7. Print a clear PASS/INFO line:
//       - if a verified best-known value exists: print achieved cut, best
//         known cut, and percentage of best-known reached
//       - if no verified value exists for this instance: print the achieved
//         cut and an explicit "(reference unavailable / unconfirmed)" notice
//         instead of fabricating or guessing a number
//
// ─────────────────────────────────────────────────────────────────────────────
// SCALING NOTE
// ─────────────────────────────────────────────────────────────────────────────
// dsb_pkg::N is fixed at compile time to 20,000 (the largest Gset instance,
// G77/G81) so that ONE build can run every Gset instance.  Smaller instances
// (e.g. 800-node G1-G67) simply use active_n=800 and leave oscillators
// 800..19999 idle and untouched (the edge list never references them, so
// their CSR rows have degree 0 and are skipped by the sparse_load_controller
// in a single short SLC_CHUNK_START -> SLC_DONE pass).
// =============================================================================

`timescale 1ns/1ps

module tb_gset_benchmark;

    import dsb_pkg::*;
    import dsb_sparse_pkg::*;
    import gset_bestknown_pkg::*;

    // =========================================================================
    // Plusarg-configurable run parameters
    // =========================================================================
    string  gset_file;
    int     gset_id;
    int     nstep_cfg;

    initial begin
        if (!$value$plusargs("gset_file=%s", gset_file))
            gset_file = "G1.txt";
        if (!$value$plusargs("gset_id=%d", gset_id))
            gset_id = 1;
        if (!$value$plusargs("nstep=%d", nstep_cfg))
            nstep_cfg = 30;
    end

    // ── Schedule constants (fixed-point, Q0.14 unless noted) ──────────────────
    localparam int A0_FP = 16384;   // 1.0
    localparam int DT_FP = 164;     // ~0.01
    localparam int C0_FP = 410;     // ~0.025 (tuned down vs the 8-node demo
                                     // since average degree is much higher
                                     // at Gset scale -- keeps the force term
                                     // from dominating the wall too quickly)

    // ── Clock / reset ─────────────────────────────────────────────────────────
    logic clk   = 1'b0;
    logic rst_n = 1'b0;
    always #5 clk = ~clk;

    // =========================================================================
    // DUT port signals
    // =========================================================================
    logic                    mm_wr_en    = 1'b0;
    logic [MAIN_ADDR_W-1:0]  mm_wr_addr  = '0;
    logic [MAIN_WORD_W-1:0]  mm_wr_data  = '0;

    logic [$clog2(N)-1:0]    active_n_in = '0;
    logic [CSR_ADDR_W-1:0]   num_edges_in = '0;

    logic [$clog2(N)-1:0]    wr_xy_idx   = '0;
    logic                    wr_xy_en    = 1'b0;
    logic signed [XY_W-1:0]  wr_x        = '0;
    logic signed [XY_W-1:0]  wr_y        = '0;

    logic                    run         = 1'b0;
    logic [STEP_W-1:0]       Nstep_in;
    logic [A_BITS-1:0]       a0_fp_in;
    logic [XY_FRAC-1:0]      dt_fp_in;
    logic [XY_FRAC-1:0]      c0_fp_in;

    assign a0_fp_in  = A_BITS'(A0_FP);
    assign dt_fp_in  = XY_FRAC'(DT_FP);
    assign c0_fp_in  = XY_FRAC'(C0_FP);

    logic [N*XY_W-1:0]       x_out;
    logic [N*XY_W-1:0]       y_out;
    logic [N-1:0]             signs_out;
    logic                    step_done;
    logic                    schedule_done;
    logic                    scan_done;

    // =========================================================================
    // DUT instantiation
    // =========================================================================
    sparse_dsb_array #(
        .N_P       (N),
        .IC_BITS_P (IC_BITS),
        .XY_W_P    (XY_W),
        .XY_FRAC_P (XY_FRAC),
        .A_BITS_P  (A_BITS),
        .STEP_W_P  (STEP_W),
        .ACCUM_W_P (ACCUM_W),
        .PROD_W_P  (PROD_W)
    ) u_dut (
        .clk           (clk),
        .rst_n         (rst_n),
        .mm_wr_en      (mm_wr_en),
        .mm_wr_addr    (mm_wr_addr),
        .mm_wr_data    (mm_wr_data),
        .active_n      (active_n_in),
        .num_edges     (num_edges_in),
        .wr_xy_idx     (wr_xy_idx),
        .wr_xy_en      (wr_xy_en),
        .wr_x          (wr_x),
        .wr_y          (wr_y),
        .run           (run),
        .Nstep         (Nstep_in),
        .a0_fp         (a0_fp_in),
        .dt_fp         (dt_fp_in),
        .c0_fp         (c0_fp_in),
        .x_out         (x_out),
        .y_out         (y_out),
        .signs_out     (signs_out),
        .step_done     (step_done),
        .schedule_done (schedule_done),
        .scan_done     (scan_done)
    );

    // =========================================================================
    // Gset file storage (parsed once at time 0, held in TB-side arrays so the
    // cut value can be re-evaluated later without re-reading the file)
    // =========================================================================
    int unsigned num_nodes;
    int unsigned num_edges;

    // Parsed edge list (1-indexed in the file; stored 0-indexed here)
    int unsigned edge_row [0:MAX_EDGES-1];
    int unsigned edge_col [0:MAX_EDGES-1];
    int signed   edge_wt  [0:MAX_EDGES-1];

    // =========================================================================
    // Task: parse a Gset file
    // =========================================================================
    task automatic parse_gset_file(input string path);
        int fd;
        int r, a, b, w;
        fd = $fopen(path, "r");
        if (fd == 0)
            $fatal(1, "Could not open Gset file: %s", path);

        r = $fscanf(fd, "%d %d", num_nodes, num_edges);
        if (r != 2)
            $fatal(1, "Failed to parse header line of %s", path);

        if (num_nodes > N)
            $fatal(1, "Instance has %0d nodes but build only supports up to N=%0d",
                   num_nodes, N);
        if (num_edges > MAX_EDGES)
            $fatal(1, "Instance has %0d edges but build only supports up to MAX_EDGES=%0d",
                   num_edges, MAX_EDGES);

        for (int e = 0; e < num_edges; e++) begin
            r = $fscanf(fd, "%d %d %d", a, b, w);
            if (r != 3)
                $fatal(1, "Failed to parse edge %0d of %s", e, path);
            // Convert 1-indexed -> 0-indexed
            edge_row[e] = a - 1;
            edge_col[e] = b - 1;
            edge_wt[e]  = w;
        end

        $fclose(fd);
        $display("  Parsed %s: %0d nodes, %0d edges", path, num_nodes, num_edges);
    endtask

    // =========================================================================
    // Task: write one word into main_memory
    // =========================================================================
    task automatic mm_write(
        input logic [MAIN_ADDR_W-1:0] addr,
        input logic [MAIN_WORD_W-1:0] data
    );
        @(posedge clk);
        mm_wr_en   <= 1'b1;
        mm_wr_addr <= addr;
        mm_wr_data <= data;
        @(posedge clk);
        mm_wr_en   <= 1'b0;
    endtask

    // ── Pack an edge word: {row, col, weight} ─────────────────────────────────
    function automatic logic [MAIN_WORD_W-1:0] pack_edge(
        input int unsigned row, col,
        input int signed   wt
    );
        automatic logic [COL_IDX_W-1:0]  r, c;
        automatic logic [IC_BITS-1:0]    w7;
        r  = COL_IDX_W'(row);
        c  = COL_IDX_W'(col);
        w7 = IC_BITS'(signed'(wt));
        return MAIN_WORD_W'({r, c, w7});
    endfunction

    // ── Load every parsed edge into main_memory ───────────────────────────────
    task automatic load_edges();
        for (int e = 0; e < num_edges; e++) begin
            mm_write(MAIN_ADDR_W'(MAIN_EDGE_BASE + e),
                     pack_edge(edge_row[e], edge_col[e], edge_wt[e]));
        end
        $display("  Loaded %0d edges into main_memory.", num_edges);
    endtask

    // ── Initialise x with small alternating perturbations (standard dSB seed)
    // ── Avoids the symmetric all-equal trap where every oscillator starts
    // ── with the exact same sign and never differentiates.
    task automatic init_x_y();
        automatic int seed = 12345;
        for (int i = 0; i < num_nodes; i++) begin
            automatic int unsigned r;
            automatic logic signed [XY_W-1:0] xv;
            r  = $urandom(seed + i) % 200;       // 0..199
            xv = XY_W'(signed'(r) - 100);          // -100..+99 -> small seed
            @(posedge clk);
            wr_xy_en  <= 1'b1;
            wr_xy_idx <= $clog2(N)'(i);
            wr_x      <= xv;
            wr_y      <= '0;
            @(posedge clk);
            wr_xy_en  <= 1'b0;
        end
        $display("  Initialised x for %0d oscillators with small random seed.", num_nodes);
    endtask

    // ── Wait for a signal with timeout ────────────────────────────────────────
    task automatic wait_high(ref logic sig, input longint timeout, input string lbl);
        longint cnt = 0;
        while (!sig) begin
            @(posedge clk);
            if (++cnt > timeout)
                $fatal(1, "[TIMEOUT] %s did not assert within %0d cycles", lbl, timeout);
        end
    endtask

    // =========================================================================
    // Max-Cut evaluator
    // -------------------------------------------------------------------------
    // cut = sum over edges (i,j,w) of  w  if sign(x_i) != sign(x_j) else 0
    // (For Gset's {-1,+1} weighted instances this correctly counts negative
    //  contributions for edges that are NOT cut, matching the standard
    //  objective: max_x sum_ij w_ij * (1 - x_i*x_j)/2 )
    // =========================================================================
    function automatic longint compute_cut();
        longint cut = 0;
        for (int e = 0; e < num_edges; e++) begin
            automatic logic si, sj;
            si = signs_out[edge_row[e]];
            sj = signs_out[edge_col[e]];
            if (si !== sj)
                cut += edge_wt[e];
        end
        return cut;
    endfunction

    // =========================================================================
    // Stimulus
    // =========================================================================
    initial begin
        $display("");
        $display("============================================================");
        $display(" Sparse dSB Gset Benchmark");
        $display(" File   : %s", gset_file);
        $display(" G-id   : G%0d", gset_id);
        $display(" Nstep  : %0d", nstep_cfg);
        $display("============================================================");

        // ── Parse the Gset file BEFORE reset, so sizes are known ──────────────
        parse_gset_file(gset_file);
        Nstep_in = STEP_W'(nstep_cfg);

        // ── Release reset ──────────────────────────────────────────────────────
        repeat(4) @(posedge clk);
        rst_n <= 1'b1;
        repeat(2) @(posedge clk);

        // =====================================================================
        // PHASE 1 — Load edge list
        // =====================================================================
        $display("");
        $display("--- Phase 1: Loading edge list into main_memory ---");
        load_edges();

        // =====================================================================
        // PHASE 2 — Initialise x/y
        // =====================================================================
        $display("");
        $display("--- Phase 2: Initialising oscillator state ---");
        init_x_y();

        // =====================================================================
        // PHASE 3 — Configure instance size and run
        // =====================================================================
        $display("");
        $display("--- Phase 3: Running sparse dSB solver ---");
        active_n_in  <= $clog2(N)'(num_nodes);
        num_edges_in <= CSR_ADDR_W'(num_edges);

        @(posedge clk);
        run <= 1'b1;

        $display("  Waiting for edge-list scan (CSR build)...");
        wait_high(scan_done, 20_000_000, "scan_done");
        $display("  Scan complete.");

        $display("  Running %0d dSB steps over %0d oscillators...", nstep_cfg, num_nodes);
        wait_high(schedule_done, 2_000_000_000, "schedule_done");
        @(posedge clk);
        run <= 1'b0;
        $display("  schedule_done asserted.");

        // =====================================================================
        // PHASE 4 — Evaluate cut and compare to best-known
        // =====================================================================
        $display("");
        $display("--- Phase 4: Evaluating result ---");

        begin
            automatic longint achieved_cut;
            automatic gset_entry_t ref_entry;
            automatic real pct;

            achieved_cut = compute_cut();
            ref_entry    = gset_lookup(gset_id);

            $display("");
            $display("============================================================");
            $display(" RESULT  —  G%0d  (%0d nodes, %0d edges)", gset_id, num_nodes, num_edges);
            $display("============================================================");
            $display("  Achieved Max-Cut (this design) : %0d", achieved_cut);

            if (ref_entry.has_reference) begin
                pct = 100.0 * real'(achieved_cut) / real'(ref_entry.best_known_cut);
                $display("  Best-known Max-Cut (reference) : %0d", ref_entry.best_known_cut);
                $display("  Percentage of best-known       : %0.2f %%", pct);
                if (achieved_cut >= ref_entry.best_known_cut)
                    $display("  >>> MATCHES OR EXCEEDS best-known reference value! <<<");
                else
                    $display("  (within %0.2f%% of best-known)", 100.0 - pct);
            end else begin
                $display("  Best-known Max-Cut (reference) : N/A");
                $display("  >>> NOTE: no confirmed/citable best-known value is");
                $display("  >>>       available for G%0d in this testbench's reference", gset_id);
                $display("  >>>       table.  Only THIS design's result is shown above.");
                $display("  >>>       (Add a verified entry to gset_bestknown_pkg.sv");
                $display("  >>>        to enable comparison for this instance.)");
            end
            $display("============================================================");
        end

        repeat(4) @(posedge clk);
        $display("");
        $display("Benchmark run complete.");
        $finish;
    end

    // ── Global watchdog (generous -- 20,000-node sweeps take many cycles) ─────
    initial begin
        #2_000_000_000;  // 2 second sim time ceiling
        $fatal(1, "GLOBAL TIMEOUT");
    end

endmodule

// =============================================================================
// tb_random_problem.sv
//
// Generic testbench for exercising the sparse dSB design against a randomly
// generated Max-Cut / Ising problem -- NOT a Gset instance, and with NO
// best-known reference value, by design.  This is a "does my design produce
// something sane and good on an arbitrary, never-seen-before problem" check,
// not a benchmark against the literature.
//
// ─────────────────────────────────────────────────────────────────────────────
// WHAT THIS TESTBENCH DOES DIFFERENTLY FROM tb_gset_benchmark.sv
// ─────────────────────────────────────────────────────────────────────────────
//  - Generates the graph itself (Erdos-Renyi style random edges, random
//    +-1 weights) instead of parsing a Gset file. No file I/O needed.
//  - Does NOT look up any best-known value -- there isn't one; the whole
//    point is testing behaviour on a problem nobody has solved before.
//  - Reports the achieved cut alongside two FREE, problem-agnostic yardsticks
//    that exist for ANY graph, known or not:
//      (a) Random-cut baseline: the EXPECTED cut value if signs were
//          assigned uniformly at random. For an unweighted-ish edge set
//          this is ~ (sum of |weights|) / 2 in expectation; we compute the
//          exact expectation analytically rather than sampling.
//      (b) A greedy local-search baseline computed in the testbench itself
//          (pure SystemVerilog, no DUT): repeatedly flip whichever single
//          node's sign increases the cut the most, until no flip helps.
//          This is a cheap, well-understood reference algorithm (local
//          search / 1-opt) that gives a meaningful "is my design doing
//          better than a simple greedy heuristic" signal without needing
//          ANY external literature value.
//  - Sweeps several random seeds / graph sizes in one run so you can see
//    how behaviour changes with problem size and density.
//
// ─────────────────────────────────────────────────────────────────────────────
// HOW TO CONFIGURE
// ─────────────────────────────────────────────────────────────────────────────
// Via plusargs (no recompilation needed):
//   +rnd_n=<nodes>        number of nodes (default 200)
//   +rnd_avgdeg=<deg>     average degree (default 12)
//   +rnd_seed=<seed>      RNG seed (default 1)
//   +rnd_signed=<0|1>     1 = weights in {-1,+1} (default), 0 = all weights +1
//   +nstep=<N>            dSB schedule steps (default 40)
//   +rnd_trials=<N>       how many independent random problems to run in
//                         sequence, each with seed+trial_index (default 1)
//
// Example:
//   vvp sparse_dsb_sim +rnd_n=2000 +rnd_avgdeg=20 +rnd_seed=7 +rnd_trials=5
// =============================================================================

`timescale 1ns/1ps

module tb_random_problem;

    import dsb_pkg::*;
    import dsb_sparse_pkg::*;

    // =========================================================================
    // Plusarg-configurable run parameters
    // =========================================================================
    int rnd_n;
    int rnd_avgdeg;
    int rnd_seed;
    int rnd_signed;
    int nstep_cfg;
    int rnd_trials;

    initial begin
        if (!$value$plusargs("rnd_n=%d",       rnd_n))       rnd_n       = 200;
        if (!$value$plusargs("rnd_avgdeg=%d",  rnd_avgdeg))  rnd_avgdeg  = 12;
        if (!$value$plusargs("rnd_seed=%d",    rnd_seed))    rnd_seed    = 1;
        if (!$value$plusargs("rnd_signed=%d",  rnd_signed))  rnd_signed  = 1;
        if (!$value$plusargs("nstep=%d",       nstep_cfg))   nstep_cfg   = 40;
        if (!$value$plusargs("rnd_trials=%d",  rnd_trials))  rnd_trials  = 1;
    end

    // ── Schedule constants ────────────────────────────────────────────────────
    localparam int A0_FP = 16384;
    localparam int DT_FP = 164;
    localparam int C0_FP = 410;

    // ── Clock / reset ─────────────────────────────────────────────────────────
    logic clk   = 1'b0;
    logic rst_n = 1'b0;
    always #5 clk = ~clk;

    // =========================================================================
    // DUT port signals (identical wiring to tb_gset_benchmark.sv)
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

    sparse_dsb_array #(
        .N_P       (N), .IC_BITS_P (IC_BITS), .XY_W_P (XY_W),
        .XY_FRAC_P (XY_FRAC), .A_BITS_P (A_BITS), .STEP_W_P (STEP_W),
        .ACCUM_W_P (ACCUM_W), .PROD_W_P (PROD_W)
    ) u_dut (
        .clk           (clk),         .rst_n         (rst_n),
        .mm_wr_en      (mm_wr_en),    .mm_wr_addr    (mm_wr_addr),
        .mm_wr_data    (mm_wr_data),
        .active_n      (active_n_in), .num_edges     (num_edges_in),
        .wr_xy_idx     (wr_xy_idx),   .wr_xy_en      (wr_xy_en),
        .wr_x          (wr_x),        .wr_y          (wr_y),
        .run           (run),         .Nstep         (Nstep_in),
        .a0_fp         (a0_fp_in),    .dt_fp         (dt_fp_in),
        .c0_fp         (c0_fp_in),
        .x_out         (x_out),       .y_out         (y_out),
        .signs_out     (signs_out),   .step_done     (step_done),
        .schedule_done (schedule_done),
        .scan_done     (scan_done)
    );

    // =========================================================================
    // Random problem storage (testbench-side, regenerated per trial)
    // =========================================================================
    int unsigned num_nodes;
    int unsigned num_edges;
    int unsigned edge_row [0:MAX_EDGES-1];
    int unsigned edge_col [0:MAX_EDGES-1];
    int signed   edge_wt  [0:MAX_EDGES-1];

    // =========================================================================
    // Task: generate a random Erdos-Renyi-style graph
    // -------------------------------------------------------------------------
    // Picks (n * avg_deg / 2) distinct undirected edges uniformly at random,
    // each with a random weight from {-1,+1} (or always +1 if rnd_signed=0).
    // =========================================================================
    task automatic generate_random_problem(
        input int n, avg_deg, seed, use_signed
    );
        automatic int target_edges;
        automatic int placed;
        automatic bit seen [int];   // associative array: dedupe edges

        num_nodes = n;
        target_edges = (n * avg_deg) / 2;
        if (target_edges > MAX_EDGES) begin
            $display("  WARNING: requested %0d edges exceeds MAX_EDGES=%0d -- clamping",
                      target_edges, MAX_EDGES);
            target_edges = MAX_EDGES;
        end

        placed = 0;
        seen.delete();
        while (placed < target_edges) begin
            automatic int a, b, key;
            automatic int lo, hi;
            a = $urandom(seed + placed*2)     % n;
            b = $urandom(seed + placed*2 + 1) % n;
            if (a == b) continue;
            lo = (a < b) ? a : b;
            hi = (a < b) ? b : a;
            key = lo * n + hi;
            if (seen.exists(key)) continue;
            seen[key] = 1;

            edge_row[placed] = lo;
            edge_col[placed] = hi;
            edge_wt[placed]  = use_signed
                                ? (($urandom(seed + 9999 + placed) % 2) ? 1 : -1)
                                : 1;
            placed++;
        end
        num_edges = placed;

        $display("  Generated random problem: %0d nodes, %0d edges, avg_deg=%0.1f, weights=%s",
                  num_nodes, num_edges, 2.0*num_edges/num_nodes,
                  use_signed ? "{-1,+1}" : "{+1}");
    endtask

    // =========================================================================
    // Tasks: load/init (identical pattern to tb_gset_benchmark.sv)
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

    task automatic load_edges();
        for (int e = 0; e < num_edges; e++)
            mm_write(MAIN_ADDR_W'(MAIN_EDGE_BASE + e),
                     pack_edge(edge_row[e], edge_col[e], edge_wt[e]));
    endtask

    task automatic init_x_y(input int seed);
        for (int i = 0; i < num_nodes; i++) begin
            automatic int unsigned r;
            automatic logic signed [XY_W-1:0] xv;
            r  = $urandom(seed + 555 + i) % 200;
            xv = XY_W'(signed'(r) - 100);
            @(posedge clk);
            wr_xy_en  <= 1'b1;
            wr_xy_idx <= $clog2(N)'(i);
            wr_x      <= xv;
            wr_y      <= '0;
            @(posedge clk);
            wr_xy_en  <= 1'b0;
        end
    endtask

    task automatic wait_high(ref logic sig, input longint timeout, input string lbl);
        longint cnt = 0;
        while (!sig) begin
            @(posedge clk);
            if (++cnt > timeout)
                $fatal(1, "[TIMEOUT] %s did not assert within %0d cycles", lbl, timeout);
        end
    endtask

    // =========================================================================
    // Max-Cut evaluator (same definition as tb_gset_benchmark.sv)
    // =========================================================================
    function automatic longint compute_cut_from_signs(ref logic [N-1:0] signs);
        longint cut = 0;
        for (int e = 0; e < num_edges; e++) begin
            if (signs[edge_row[e]] !== signs[edge_col[e]])
                cut += edge_wt[e];
        end
        return cut;
    endfunction

    // =========================================================================
    // FREE BASELINE 1 — analytical random-cut expectation
    // -------------------------------------------------------------------------
    // For independent uniform-random +-1 sign assignment, each edge is cut
    // with probability 1/2 regardless of its weight, so:
    //   E[cut] = sum_e ( w_e * 0.5 )   [edges contribute w_e when cut, 0 when not]
    // This is an exact closed-form expectation, not a sampled estimate, and
    // exists for ANY graph -- known or brand new -- making it a always-available
    // sanity floor: a competent solver should comfortably beat this.
    // =========================================================================
    function automatic real expected_random_cut();
        real sum_w = 0.0;
        for (int e = 0; e < num_edges; e++)
            sum_w += real'(edge_wt[e]);
        return 0.5 * sum_w;
    endfunction

    // =========================================================================
    // FREE BASELINE 2 — greedy 1-opt local search (computed in the TB itself)
    // -------------------------------------------------------------------------
    // Standard, well-understood baseline: start from a random sign vector,
    // repeatedly flip whichever single node gives the largest cut increase,
    // stop when no flip helps (local optimum). This needs no external
    // reference and works for any graph. Implemented with an adjacency list
    // built on the fly (associative array of dynamic arrays) purely in the
    // testbench -- this code never touches the DUT.
    // =========================================================================
    function automatic longint greedy_local_search(input int seed);
        // Build adjacency: for each node, list of (neighbour, weight)
        typedef struct { int neighbour; int weight; } adj_t;
        adj_t adj_list [int][$];   // associative array of queues, keyed by node

        bit sign_arr [int];        // 1 = positive, 0 = negative
        longint cut;
        bit improved;

        adj_list.delete();
        for (int e = 0; e < num_edges; e++) begin
            adj_list[edge_row[e]].push_back('{edge_col[e], edge_wt[e]});
            adj_list[edge_col[e]].push_back('{edge_row[e], edge_wt[e]});
        end

        // Random initial assignment
        for (int i = 0; i < num_nodes; i++)
            sign_arr[i] = $urandom(seed + 777 + i) % 2;

        // Iterate to local optimum (bounded iterations to guarantee termination
        // in simulation time even on pathological cases)
        for (int pass = 0; pass < 50; pass++) begin
            improved = 0;
            for (int i = 0; i < num_nodes; i++) begin
                automatic longint delta = 0;
                // Flipping node i changes the cut status of every incident edge.
                // delta = (new cut contribution) - (old cut contribution)
                //       = sum over neighbours of: w * (+1 if now-cut else -1)
                //         relative to current state
                foreach (adj_list[i][k]) begin
                    automatic bit same_before;
                    same_before = (sign_arr[i] == sign_arr[adj_list[i][k].neighbour]);
                    // Before flip: cut contributes w if !same_before else 0
                    // After flip : same_before inverts, so contribution inverts
                    if (same_before)
                        delta += adj_list[i][k].weight;   // was not cut, becomes cut
                    else
                        delta -= adj_list[i][k].weight;   // was cut, becomes not cut
                end
                if (delta > 0) begin
                    sign_arr[i] = ~sign_arr[i];
                    improved = 1;
                end
            end
            if (!improved) break;
        end

        // Evaluate final cut
        cut = 0;
        for (int e = 0; e < num_edges; e++)
            if (sign_arr[edge_row[e]] !== sign_arr[edge_col[e]])
                cut += edge_wt[e];

        return cut;
    endfunction

    // =========================================================================
    // Run one full trial: generate -> load -> run DUT -> evaluate -> report
    // =========================================================================
    task automatic run_trial(input int trial_idx);
        automatic int seed_this_trial = rnd_seed + trial_idx * 1000;
        automatic longint dut_cut, greedy_cut;
        automatic real rand_cut, pct_vs_greedy;

        $display("");
        $display("============================================================");
        $display(" Trial %0d  (seed=%0d)", trial_idx, seed_this_trial);
        $display("============================================================");

        generate_random_problem(rnd_n, rnd_avgdeg, seed_this_trial, rnd_signed);

        // ── Reset DUT state for this trial ────────────────────────────────────
        rst_n <= 1'b0;
        repeat(4) @(posedge clk);
        rst_n <= 1'b1;
        repeat(2) @(posedge clk);

        $display("  Loading edges into main_memory...");
        load_edges();

        $display("  Initialising oscillator state...");
        init_x_y(seed_this_trial);

        active_n_in  <= $clog2(N)'(num_nodes);
        num_edges_in <= CSR_ADDR_W'(num_edges);
        Nstep_in     <= STEP_W'(nstep_cfg);

        @(posedge clk);
        run <= 1'b1;

        wait_high(scan_done, 20_000_000, "scan_done");
        wait_high(schedule_done, 2_000_000_000, "schedule_done");
        @(posedge clk);
        run <= 1'b0;

        dut_cut    = compute_cut_from_signs(signs_out);
        rand_cut   = expected_random_cut();
        greedy_cut = greedy_local_search(seed_this_trial);

        $display("");
        $display("  --- Result ---");
        $display("  Sparse dSB design cut value      : %0d", dut_cut);
        $display("  Random-assignment expected cut   : %0.1f  (analytical floor)", rand_cut);
        $display("  Greedy 1-opt local search cut     : %0d  (cheap reference algorithm)", greedy_cut);

        if (real'(dut_cut) > rand_cut)
            $display("  PASS  design beats the random-assignment floor (%0d > %0.1f)", dut_cut, rand_cut);
        else
            $display("  FAIL  design did NOT beat the random-assignment floor (%0d <= %0.1f)", dut_cut, rand_cut);

        if (greedy_cut != 0)
            pct_vs_greedy = 100.0 * real'(dut_cut) / real'(greedy_cut);
        else
            pct_vs_greedy = 0.0;

        if (dut_cut >= greedy_cut)
            $display("  PASS  design matches or beats greedy local search (%0.1f%% of greedy)", pct_vs_greedy);
        else
            $display("  INFO  design is below greedy local search (%0.1f%% of greedy) -- consider", pct_vs_greedy);
        $display("        more dSB steps (current Nstep=%0d) or schedule tuning", nstep_cfg);
    endtask

    // =========================================================================
    // Stimulus — run rnd_trials independent random problems
    // =========================================================================
    initial begin
        $display("");
        $display("============================================================");
        $display(" Sparse dSB Random-Problem Stress Test");
        $display(" (no best-known reference -- evaluating against problem-");
        $display("  agnostic baselines only, since the problem is novel)");
        $display("============================================================");
        $display(" nodes     : %0d", rnd_n);
        $display(" avg degree: %0d", rnd_avgdeg);
        $display(" weights   : %s", rnd_signed ? "{-1,+1}" : "{+1}");
        $display(" Nstep     : %0d", nstep_cfg);
        $display(" trials    : %0d", rnd_trials);

        for (int t = 0; t < rnd_trials; t++)
            run_trial(t);

        $display("");
        $display("============================================================");
        $display(" All %0d trial(s) complete.", rnd_trials);
        $display("============================================================");
        $finish;
    end

    initial begin
        #2_000_000_000;
        $fatal(1, "GLOBAL TIMEOUT");
    end

endmodule

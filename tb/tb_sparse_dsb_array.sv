// tb_sparse_dsb_array.sv — multi‑case test with correct declarations
`timescale 1ns/1ps

module tb_sparse_dsb_array;

    import dsb_pkg::*;
    import dsb_sparse_pkg::*;

    // ── Clock / reset ──────────────────────────────────────────────────────
    logic clk = 0;
    always #5 clk = ~clk;

    logic rst_n = 0;

    // ── Top‑level DUT signals ─────────────────────────────────────────────
    logic                    mm_wr_en;
    logic [MAIN_ADDR_W-1:0]  mm_wr_addr;
    logic [MAIN_WORD_W-1:0]  mm_wr_data;

    logic [$clog2(N)-1:0]    active_n;
    logic [CSR_ADDR_W-1:0]   num_edges;

    logic [$clog2(N)-1:0]    wr_xy_idx;
    logic                    wr_xy_en;
    logic signed [XY_W-1:0]  wr_x;
    logic signed [XY_W-1:0]  wr_y;

    logic                    run;
    logic [STEP_W-1:0]       Nstep;
    logic [A_BITS-1:0]       a0_fp;
    logic [XY_FRAC-1:0]      dt_fp;
    logic [XY_FRAC-1:0]      c0_fp;

    logic [N*XY_W-1:0]       x_out;
    logic [N*XY_W-1:0]       y_out;
    logic [N-1:0]            signs_out;
    logic                    step_done;
    logic                    schedule_done;
    logic                    scan_done;

    // DUT
    sparse_dsb_array #(
        .N_P        (N),
        .IC_BITS_P  (IC_BITS),
        .XY_W_P     (XY_W),
        .XY_FRAC_P  (XY_FRAC),
        .A_BITS_P   (A_BITS),
        .STEP_W_P   (STEP_W),
        .ACCUM_W_P  (ACCUM_W),
        .PROD_W_P   (PROD_W)
    ) u_dut (
        .clk         (clk),
        .rst_n       (rst_n),
        .mm_wr_en    (mm_wr_en),
        .mm_wr_addr  (mm_wr_addr),
        .mm_wr_data  (mm_wr_data),
        .active_n    (active_n),
        .num_edges   (num_edges),
        .wr_xy_idx   (wr_xy_idx),
        .wr_xy_en    (wr_xy_en),
        .wr_x        (wr_x),
        .wr_y        (wr_y),
        .run         (run),
        .Nstep       (Nstep),
        .a0_fp       (a0_fp),
        .dt_fp       (dt_fp),
        .c0_fp       (c0_fp),
        .x_out       (x_out),
        .y_out       (y_out),
        .signs_out   (signs_out),
        .step_done   (step_done),
        .schedule_done(schedule_done),
        .scan_done   (scan_done)
    );

    // ── Constants ─────────────────────────────────────────────────────────
    localparam int A0_FP = 16384;
    localparam int DT_FP = 164;
    localparam int C0_FP = 8192;

    // ── Testbench helpers ──────────────────────────────────────────────────
    int step_count = 0;

    task automatic mm_write(input logic [MAIN_ADDR_W-1:0] addr,
                            input logic [MAIN_WORD_W-1:0] data);
        @(posedge clk);
        mm_wr_en = 1;
        mm_wr_addr = addr;
        mm_wr_data = data;
        @(posedge clk);
        mm_wr_en = 0;
    endtask

    function automatic logic [MAIN_WORD_W-1:0] pack_edge(input int r, c, w);
        return {COL_IDX_W'(r), COL_IDX_W'(c), IC_BITS'(signed'(w))};
    endfunction

    task automatic write_x(input int idx, input logic signed [XY_W-1:0] val);
        @(posedge clk);
        wr_xy_en = 1;
        wr_xy_idx = idx;
        wr_x = val;
        wr_y = '0;
        @(posedge clk);
        wr_xy_en = 0;
    endtask

    task automatic wait_high(ref logic sig, input int timeout, input string lbl);
        int cnt = 0;
        while (!sig) begin
            @(posedge clk);
            cnt++;
            if (cnt > timeout) $fatal(1, "TIMEOUT: %s", lbl);
        end
    endtask

    // ── Test case structure ───────────────────────────────────────────────
    typedef struct {
        string   name;
        int      num_nodes;
        int      num_edges;
        int      nstep;
        int      edges[ ][3];   // {src, dst, weight}
        int      x_init[ ];
        int      min_cut;       // minimum acceptable cut for this case
        int      max_energy;    // energy must be <= this
    } test_case_t;

    // ── Helper to run one test case ──────────────────────────────────────
    task automatic run_case(input test_case_t tc, output bit pass);
        int cut;
        int energy;
        int e;

        $display("");
        $display("====================================================");
        $display(" Test case: %s", tc.name);
        $display("  nodes = %0d, edges = %0d, Nstep = %0d",
                 tc.num_nodes, tc.num_edges, tc.nstep);
        $display("====================================================");

        // 1. Load edges
        for (e = 0; e < tc.num_edges; e++) begin
            int r = tc.edges[e][0];
            int c = tc.edges[e][1];
            int w = tc.edges[e][2];
            mm_write(MAIN_EDGE_BASE + e, pack_edge(r, c, w));
        end

        // 2. Initialize x
        for (int i = 0; i < tc.num_nodes; i++) begin
            write_x(i, XY_W'(tc.x_init[i] << (XY_FRAC - 4)));
        end

        // 3. Set parameters and start
        active_n  = $clog2(N)'(tc.num_nodes);
        num_edges = CSR_ADDR_W'(tc.num_edges);
        Nstep     = STEP_W'(tc.nstep);
        a0_fp     = A_BITS'(A0_FP);
        dt_fp     = XY_FRAC'(DT_FP);
        c0_fp     = XY_FRAC'(C0_FP);

        @(posedge clk);
        run = 1;

        wait_high(scan_done, 10000, "scan_done");
        wait_high(schedule_done, 1_000_000, "schedule_done");
        @(posedge clk);
        run = 0;

        // 4. Evaluate result
        cut = 0;
        for (e = 0; e < tc.num_edges; e++) begin
            int r = tc.edges[e][0];
            int c = tc.edges[e][1];
            if (signs_out[r] !== signs_out[c]) cut++;
        end

        // Compute energy
        energy = 0;
        for (e = 0; e < tc.num_edges; e++) begin
            int r = tc.edges[e][0];
            int c = tc.edges[e][1];
            int w = tc.edges[e][2];
            int s1 = signs_out[r] ? 1 : -1;
            int s2 = signs_out[c] ? 1 : -1;
            energy += w * s1 * s2;
        end

        $display("  Achieved cut = %0d / %0d", cut, tc.num_edges);
        $display("  Energy = %0d", energy);

        pass = 1'b1;
        if (cut < tc.min_cut) begin
            $display("  FAIL: cut %0d < min required %0d", cut, tc.min_cut);
            pass = 1'b0;
        end
        if (energy > tc.max_energy) begin
            $display("  FAIL: energy %0d > max allowed %0d", energy, tc.max_energy);
            pass = 1'b0;
        end
        if (cut == tc.num_edges) begin
            $display("  OPTIMAL: all edges cut!");
        end

        if (pass) $display("  PASS");
        else      $display("  FAIL");
        $display("====================================================");
        #100;
    endtask

    // ── Test case generators ─────────────────────────────────────────────
    function automatic test_case_t cycle_case(int n, int nstep);
        test_case_t tc;
        tc.name = $sformatf("%0d‑node cycle", n);
        tc.num_nodes = n;
        tc.num_edges = n;
        tc.nstep = nstep;
        tc.edges = new[n];
        for (int i = 0; i < n; i++) begin
            int j = (i + 1) % n;
            tc.edges[i] = '{i, j, -1};
        end
        tc.x_init = new[n];
        for (int i = 0; i < n; i++) tc.x_init[i] = (i % 2 == 0) ? 50 : -30;
        tc.min_cut = n/2;
        tc.max_energy = -1;
        return tc;
    endfunction

    function automatic test_case_t complete_bipartite_case(int a, int b, int nstep);
        test_case_t tc;
        int n = a + b;
        int idx = 0;
        tc.name = $sformatf("K%0d,%0d complete bipartite", a, b);
        tc.num_nodes = n;
        tc.num_edges = a * b;
        tc.nstep = nstep;
        tc.edges = new[tc.num_edges];
        for (int i = 0; i < a; i++) begin
            for (int j = a; j < n; j++) begin
                tc.edges[idx++] = '{i, j, -1};
            end
        end
        tc.x_init = new[n];
        for (int i = 0; i < n; i++) tc.x_init[i] = (i < a) ? 40 : -40;
        tc.min_cut = tc.num_edges / 2;
        tc.max_energy = -1;
        return tc;
    endfunction

    function automatic test_case_t random_graph_case(int n, int density, int nstep);
        test_case_t tc;
        int max_edges, num_edges, idx;
        tc.name = $sformatf("random graph n=%0d density=%0.1f", n, density/10.0);
        tc.num_nodes = n;
        max_edges = n*(n-1)/2;
        num_edges = max_edges * density / 10;
        tc.num_edges = num_edges;
        tc.nstep = nstep;
        tc.edges = new[num_edges];
        idx = 0;
        for (int i = 0; i < n && idx < num_edges; i++) begin
            for (int j = i+1; j < n && idx < num_edges; j++) begin
                if ( (i*j + j + i) % 10 < density ) begin
                    tc.edges[idx++] = '{i, j, -1};
                end
            end
        end
        // Pad if necessary
        while (idx < num_edges) begin
            int i = idx % n;
            int j = (idx + 7) % n;
            if (i != j) begin
                tc.edges[idx++] = '{i, j, -1};
            end
        end
        tc.x_init = new[n];
        for (int i = 0; i < n; i++) tc.x_init[i] = (i % 2 == 0) ? 30 : -30;
        tc.min_cut = num_edges / 3;
        tc.max_energy = -1;
        return tc;
    endfunction

    // ── Main test sequence ──────────────────────────────────────────────
    initial begin
        // ---- All declarations at the top ----
        test_case_t cases[4];
        int total_passed;
        int total_cases;
        bit pass;

        $dumpfile("tb_sparse_dsb_array.vcd");
        $dumpvars(0, tb_sparse_dsb_array);

        // Build test case list
        cases[0] = cycle_case(4, 200);
        cases[1] = cycle_case(8, 300);
        cases[2] = complete_bipartite_case(4,4, 400);
        cases[3] = random_graph_case(10, 3, 500);

        rst_n = 0;
        repeat (4) @(posedge clk);
        rst_n = 1;
        repeat (2) @(posedge clk);
        $display("[%0t] Reset released", $time);

        total_passed = 0;
        total_cases = $size(cases);   // <-- FIXED: use $size() for static array

        for (int c = 0; c < total_cases; c++) begin
            run_case(cases[c], pass);
            if (pass) total_passed++;
        end

        $display("");
        $display("====================================================");
        $display(" SUMMARY: %0d / %0d test cases PASSED", total_passed, total_cases);
        $display("====================================================");
        if (total_passed == total_cases)
            $display("ALL TESTS PASSED");
        else
            $display("SOME TESTS FAILED");
        $finish;
    end

    // ── Watchdog timeout ──────────────────────────────────────────────────
    initial begin
        #50_000_000;
        $fatal(1, "GLOBAL TIMEOUT");
    end

    // ── Step monitor ────────────────────────────────────────────────────
    always @(posedge clk) begin
        if (step_done) begin
            step_count++;
            $display("[%0t] step %0d done  signs = %0b", $time, step_count, signs_out);
        end
    end

endmodule
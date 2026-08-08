// tb_sparse_dsb_array_100nodes.sv
// Full‑chip test on K_50,50 (100 nodes, 2500 edges, optimal cut = 2500)
`timescale 1ns/1ps

module tb_sparse_dsb_array_100nodes;

    import dsb_pkg::*;
    import dsb_sparse_pkg::*;

    // ── Clock / reset ──────────────────────────────────────────────────────
    logic clk = 0;
    always #5 clk = ~clk;

    logic rst_n = 0;

    // ── DUT signals ──────────────────────────────────────────────────────
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

    // ── DUT ──────────────────────────────────────────────────────────────
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

    // ── Problem definition ───────────────────────────────────────────────
    localparam int NODES   = 100;
    localparam int HALF    = NODES / 2;   // 50
    localparam int EDGES   = HALF * HALF; // 2500
    localparam int NSTEP   = 6000;

    // Parameters for dynamics (aggressive)
    localparam int A0_FP  = 0;
    localparam int C0_FP  = 20 * 16384;   // 20.0 in Q0.14
    localparam int DT_FP  = 1638;         // 0.1 in Q0.14

    // ── Tasks ─────────────────────────────────────────────────────────────
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

    // ── Main test ────────────────────────────────────────────────────────
    initial begin
        // ---- All declarations at the top ----
        int cut, energy, e, i, r, c, s1, s2, val, edge_idx;
        bit pass;

        $dumpfile("tb_sparse_dsb_array_100nodes.vcd");
        $dumpvars(0, tb_sparse_dsb_array_100nodes);

        $display("==========================================================");
        $display(" Test: K%0d,%0d complete bipartite graph", HALF, HALF);
        $display(" Nodes = %0d, Edges = %0d, Optimal cut = %0d",
                 NODES, EDGES, EDGES);
        $display("==========================================================");

        // Reset
        rst_n = 0;
        repeat (4) @(posedge clk);
        rst_n = 1;
        repeat (2) @(posedge clk);
        $display("[%0t] Reset released", $time);

        // ── 1. Load all edges (J = +1) ──────────────────────────────────
        $display("[%0t] Loading %0d edges...", $time, EDGES);
        edge_idx = 0;
        while (edge_idx < EDGES) begin
            // Generate K_{50,50} edges: from nodes 0..49 to nodes 50..99
            r = edge_idx / HALF;        // 0..49
            c = HALF + (edge_idx % HALF); // 50..99
            mm_write(MAIN_EDGE_BASE + edge_idx, pack_edge(r, c, 1));
            if (edge_idx % 500 == 0) $display("  loaded %0d edges", edge_idx);
            edge_idx++;
        end
        $display("[%0t] All edges loaded", $time);

        // ── 2. Initialize x with small random values ──────────────────
        $display("[%0t] Initialising x with small random values", $time);
        for (i = 0; i < NODES; i++) begin
            val = (i % 2 == 0) ? 5 : -5;
            write_x(i, XY_W'(val << (XY_FRAC - 4)));
        end

        // ── 3. Set parameters and start ───────────────────────────────
        active_n  = $clog2(N)'(NODES);
        num_edges = CSR_ADDR_W'(EDGES);
        Nstep     = STEP_W'(NSTEP);
        a0_fp     = A_BITS'(A0_FP);
        dt_fp     = XY_FRAC'(DT_FP);
        c0_fp     = XY_FRAC'(C0_FP);

        $display("[%0t] Starting run with Nstep = %0d", $time, NSTEP);
        run = 1;

        wait_high(scan_done, 10000, "scan_done");
        $display("[%0t] scan_done", $time);

        wait_high(schedule_done, 20_000_000, "schedule_done");
        @(posedge clk);
        run = 0;
        $display("[%0t] schedule_done", $time);

        // ── 4. Evaluate result ────────────────────────────────────────
        cut = 0;
        for (e = 0; e < EDGES; e++) begin
            r = e / HALF;
            c = HALF + (e % HALF);
            if (signs_out[r] !== signs_out[c]) cut++;
        end

        energy = 0;
        for (e = 0; e < EDGES; e++) begin
            r = e / HALF;
            c = HALF + (e % HALF);
            s1 = signs_out[r] ? 1 : -1;
            s2 = signs_out[c] ? 1 : -1;
            energy += s1 * s2;   // J = +1
        end

        $display("==========================================================");
        $display(" FINAL STATE");
        $display("  Achieved cut = %0d / %0d", cut, EDGES);
        $display("  Energy = %0d (optimal = -%0d)", energy, EDGES);
        $display("==========================================================");

        pass = 1'b1;
        if (cut < EDGES) begin
            $display("  FAIL: cut = %0d, expected %0d", cut, EDGES);
            pass = 1'b0;
        end
        if (energy != -EDGES) begin
            $display("  FAIL: energy = %0d, expected -%0d", energy, EDGES);
            pass = 1'b0;
        end
        if (cut == EDGES && energy == -EDGES) begin
            $display("  PASS: Optimal solution found!");
        end

        // ── 5. Summary ──────────────────────────────────────────────────
        if (pass)
            $display("TEST PASSED");
        else
            $display("TEST FAILED");
        $finish;
    end

    // ── Watchdog timeout ──────────────────────────────────────────────────
    initial begin
        #2_000_000_000;
        $fatal(1, "GLOBAL TIMEOUT");
    end

    // ── Step monitor ────────────────────────────────────────────────────
    int step_count = 0;
    always @(posedge clk) begin
        if (step_done) begin
            step_count++;
            if (step_count % 100 == 0)
                $display("[%0t] step %0d done", $time, step_count);
        end
    end

endmodule
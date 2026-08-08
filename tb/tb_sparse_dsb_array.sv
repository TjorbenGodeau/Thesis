// tb_sparse_dsb_array.sv
// Self-checking testbench for sparse_dsb_array top-level module
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

    // ── Parameters for the test problem ──────────────────────────────────
    localparam int NUM_NODES = 4;
    localparam int NUM_EDGES = 4;
    localparam int NSTEP     = 200;

    localparam int A0_FP  = 16384;   // 1.0 in Q0.14
    localparam int DT_FP  = 164;     // 0.01
    localparam int C0_FP  = 8192;    // 0.5

    typedef struct { int r, c, w; } edge_t;
    edge_t EDGES [0:NUM_EDGES-1] = '{
        '{0,1,-1},
        '{1,2,-1},
        '{2,3,-1},
        '{3,0,-1}
    };

    logic signed [XY_W-1:0] X_INIT [0:NUM_NODES-1] = '{
        16'sd50,    // node0: +50
        -16'sd30,   // node1: -30
        16'sd40,    // node2: +40
        -16'sd35    // node3: -35
    };

    // ── Module‑level variable for step monitoring ──────────────────────
    int step_count = 0;

    // ── Testbench tasks ──────────────────────────────────────────────────
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

    // ── Compute Hamiltonian energy from signs ────────────────────────────
    function automatic int compute_energy(input logic [NUM_NODES-1:0] signs);
        int energy = 0;
        for (int e = 0; e < NUM_EDGES; e++) begin
            int s1 = signs[EDGES[e].r] ? 1 : -1;
            int s2 = signs[EDGES[e].c] ? 1 : -1;
            energy += EDGES[e].w * s1 * s2;
        end
        return energy;
    endfunction

    // ── Main test sequence ──────────────────────────────────────────────
    initial begin
        // ---- All declarations at the top ----
        int cut;
        int energy;
        logic success;

        $dumpfile("tb_sparse_dsb_array.vcd");
        $dumpvars(0, tb_sparse_dsb_array);

        $display("====================================================");
        $display(" Test: sparse_dsb_array on a 4‑node cycle");
        $display(" Optimal cut = %0d (all edges cut)", NUM_EDGES);
        $display("====================================================");

        rst_n = 0;
        repeat (4) @(posedge clk);
        rst_n = 1;
        repeat (2) @(posedge clk);
        $display("[%0t] Reset released", $time);

        // ── Phase 1: Load edges ───────────────────────────────────────────
        $display("[%0t] Loading %0d edges", $time, NUM_EDGES);
        for (int e = 0; e < NUM_EDGES; e++) begin
            mm_write(MAIN_EDGE_BASE + e,
                     pack_edge(EDGES[e].r, EDGES[e].c, EDGES[e].w));
            $display("  edge[%0d]: %0d -- %0d  w = %0d",
                     e, EDGES[e].r, EDGES[e].c, EDGES[e].w);
        end

        // ── Phase 2: Initialise x ─────────────────────────────────────────
        $display("[%0t] Initialising x", $time);
        for (int i = 0; i < NUM_NODES; i++) begin
            write_x(i, X_INIT[i]);
            $display("  x[%0d] = %0d", i, $signed(X_INIT[i]));
        end

        // ── Phase 3: Set problem size and run ─────────────────────────────
        active_n  = $clog2(N)'(NUM_NODES);
        num_edges = CSR_ADDR_W'(NUM_EDGES);
        Nstep     = STEP_W'(NSTEP);
        a0_fp     = A_BITS'(A0_FP);
        dt_fp     = XY_FRAC'(DT_FP);
        c0_fp     = XY_FRAC'(C0_FP);

        $display("[%0t] Starting run with Nstep = %0d", $time, NSTEP);
        run = 1;

        wait_high(scan_done, 10000, "scan_done");
        $display("[%0t] scan_done", $time);

        wait_high(schedule_done, 1_000_000, "schedule_done");
        @(posedge clk);
        run = 0;
        $display("[%0t] schedule_done", $time);

        // ── Phase 4: Check results ────────────────────────────────────────
        $display("====================================================");
        $display(" FINAL STATE");
        $display("====================================================");
        $display("signs_out = %0b", signs_out);
        $display("");

        // Compute cut and energy
        cut = 0;
        for (int e = 0; e < NUM_EDGES; e++) begin
            if (signs_out[EDGES[e].r] !== signs_out[EDGES[e].c]) begin
                cut++;
                $display("  edge %0d--%0d: CUT", EDGES[e].r, EDGES[e].c);
            end else begin
                $display("  edge %0d--%0d: not cut (both %s)",
                         EDGES[e].r, EDGES[e].c,
                         signs_out[EDGES[e].r] ? "A" : "B");
            end
        end

        energy = compute_energy(signs_out);
        $display("");
        $display("  Achieved cut = %0d / %0d (optimal = %0d)", cut, NUM_EDGES, NUM_EDGES);
        $display("  Hamiltonian energy = %0d (optimal = -%0d)", energy, NUM_EDGES);

        // Check if result is reasonable
        success = 1'b1;
        if (cut < 2) begin
            $display("  FAIL: cut too low (%0d)", cut);
            success = 1'b0;
        end else if (energy > 0) begin
            $display("  FAIL: energy positive (%0d)", energy);
            success = 1'b0;
        end else begin
            $display("  PASS: result is valid (cut >= 2, energy <= 0)");
        end
        if (cut == NUM_EDGES) begin
            $display("  OPTIMAL: all edges cut!");
        end else begin
            $display("  NOTE: suboptimal (local minimum) but acceptable for deterministic solver");
        end

        $display("====================================================");
        if (success)
            $display("TEST PASSED");
        else
            $display("TEST FAILED");
        $finish;
    end

    // ── Watchdog timeout ──────────────────────────────────────────────────
    initial begin
        #10_000_000;
        $fatal(1, "GLOBAL TIMEOUT");
    end

    // ── Monitor step progress ────────────────────────────────────────────
    always @(posedge clk) begin
        if (step_done) begin
            step_count++;
            $display("[%0t] step %0d done  signs = %0b", $time, step_count, signs_out);
        end
    end

endmodule
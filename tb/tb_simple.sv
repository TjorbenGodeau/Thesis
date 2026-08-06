// =============================================================================
// tb_simple.sv
//
// Minimal testbench for waveform inspection on a hand-verifiable 8-node
// MaxCut problem.  Designed so every signal can be traced manually.
//
// THE PROBLEM
// -----------
// 8 nodes, 10 edges, all weights +1.  This is a bipartite graph:
//
//         0
//        /|\
//       1 2 7
//       | |  \
//       3-+   6
//       |    / \
//       4   5   (back to 3,7)
//
// More precisely (0-indexed):
//   0--1, 0--2, 0--7
//   1--3, 2--3
//   3--4, 3--5
//   4--6, 5--6
//   6--7
//
// KNOWN OPTIMAL SOLUTION
// ----------------------
// Partition A = {0, 3, 6}  (sign = +1)
// Partition B = {1, 2, 4, 5, 7}  (sign = -1)
// Every single edge crosses the partition -> optimal cut = 10 (all edges)
//
// WHAT TO CHECK IN THE WAVEFORM  (see checklist at the bottom of this file)
// -------------------------------------------------------------------------
// Phase 1: mm_wr_en pulses 10 times, one per edge (20 clocks total)
// Phase 2: wr_xy_en pulses 8 times, one per oscillator
// Phase 3: scan_done goes high after ~80 clocks (8 nodes * overhead + 10 edges * 7)
// Phase 4: for each oscillator update you should see:
//            slc_start pulse -> slc_done pulse some cycles later
//            ph2_start pulse -> ph2_done pulse one cycle later
//            update_unit start -> done one cycle later
// Phase 5: schedule_done goes high after Nstep sweeps
// Phase 6: signs_out[0,3,6] = 1  signs_out[1,2,4,5,7] = 0  (optimal)
//          cut value printed = 10
//
// HAND-VERIFIABLE Jx_i VALUES (at the optimal sign vector)
// ---------------------------------------------------------
//   Jx_0 = sign(x1)+sign(x2)+sign(x7) = -1-1-1 = -3
//   Jx_1 = sign(x0)+sign(x3)          = +1+1   = +2
//   Jx_2 = sign(x0)+sign(x3)          = +1+1   = +2
//   Jx_3 = sign(x1)+sign(x2)+sign(x4)+sign(x5) = -1-1-1-1 = -4
//   Jx_4 = sign(x3)+sign(x6)          = +1+1   = +2
//   Jx_5 = sign(x3)+sign(x6)          = +1+1   = +2
//   Jx_6 = sign(x4)+sign(x5)+sign(x7) = -1-1-1 = -3
//   Jx_7 = sign(x6)+sign(x0)          = +1+1   = +2  (wait -- x0=+ and x6=+)
//          actually Jx_7 = sign(x6)+sign(x0) = +1+1 = +2... but sign_7=-1
//          force = c0*Jx_7 = positive, pushing y_7 positive -> x_7 goes positive
//          WAIT: for the optimal solution sign_7 should be -1 (partition B)
//          The force c0*Jx_7=+2 pushes x_7 TOWARD +1, but we want x_7 negative.
//          This means node 7 is a contested node -- the algorithm may not
//          always land on the global optimum for this node. That is expected
//          behaviour for dSB (it is a heuristic, not an exact solver).
//          In practice the pump term a(m)*y_7 competes with the coupling force
//          and the final sign depends on which wins. Watch x_7 in the waveform.
//
// PARAMETER CHOICES FOR THIS RUN
// --------------------------------
//   N      = 8 in dsb_pkg  (YOU MUST SET THIS -- see note below)
//   Nstep  = 20  (very short -- watch convergence)
//   dt_fp  = 164 (~0.01 in Q0.14)
//   c0_fp  = 2048 (= 0.5/avg_deg * 2^14, avg_deg=2.5 -> c0=0.2 -> fp=3277;
//                  but we use 2048 = ~0.125 to keep forces gentle for inspection)
//   a0_fp  = 16384 (= 1.0 in Q0.14)
//
// NOTE ON N
// ---------
// dsb_pkg::N is currently 20000.  For a clean short waveform, recompile with
// N=8 by temporarily changing the parameter in allModules.sv.  This shrinks
// all arrays from 20000 to 8 entries and makes the waveform readable.
// After testing, restore N=20000 for real Gset runs.
// =============================================================================

`timescale 1ns/1ps

module tb_simple;

    import dsb_pkg::*;
    import dsb_sparse_pkg::*;

    // ── Clock / reset ─────────────────────────────────────────────────────────
    logic clk   = 1'b0;
    logic rst_n = 1'b0;
    always #5 clk = ~clk;    // 100 MHz, 10 ns period

    // ── Schedule parameters ───────────────────────────────────────────────────
    // avg_degree = 2*10/8 = 2.5
    // c0 = alpha/avg_deg = 0.5/2.5 = 0.2
    // c0_fp = round(0.2 * 2^14) = round(3276.8) = 3277
    localparam int NSTEP  = 20;
    localparam int A0_FP  = 16384;   // 1.0 in Q0.14
    localparam int DT_FP  = 164;     // 0.01 in Q0.14
    localparam int C0_FP  = 3277;    // 0.2 in Q0.14

    // ── Problem size ──────────────────────────────────────────────────────────
    localparam int NUM_NODES = 8;
    localparam int NUM_EDGES = 10;

    // ── DUT signals ───────────────────────────────────────────────────────────
    logic                    mm_wr_en    = 1'b0;
    logic [MAIN_ADDR_W-1:0]  mm_wr_addr  = '0;
    logic [MAIN_WORD_W-1:0]  mm_wr_data  = '0;

    logic [$clog2(N)-1:0]    active_n_in  = '0;
    logic [CSR_ADDR_W-1:0]   num_edges_in = '0;

    logic [$clog2(N)-1:0]    wr_xy_idx = '0;
    logic                    wr_xy_en  = 1'b0;
    logic signed [XY_W-1:0]  wr_x      = '0;
    logic signed [XY_W-1:0]  wr_y      = '0;

    logic                    run      = 1'b0;
    logic [STEP_W-1:0]       Nstep_in;
    logic [A_BITS-1:0]       a0_fp_in;
    logic [XY_FRAC-1:0]      dt_fp_in;
    logic [XY_FRAC-1:0]      c0_fp_in;

    assign Nstep_in  = STEP_W'(NSTEP);
    assign a0_fp_in  = A_BITS'(A0_FP);
    assign dt_fp_in  = XY_FRAC'(DT_FP);
    assign c0_fp_in  = XY_FRAC'(C0_FP);

    logic [N*XY_W-1:0]  x_out;
    logic [N*XY_W-1:0]  y_out;
    logic [N-1:0]        signs_out;
    logic                step_done;
    logic                schedule_done;
    logic                scan_done;

    // ── DUT ───────────────────────────────────────────────────────────────────
    sparse_dsb_array #(
        .N_P(N), .IC_BITS_P(IC_BITS), .XY_W_P(XY_W), .XY_FRAC_P(XY_FRAC),
        .A_BITS_P(A_BITS), .STEP_W_P(STEP_W), .ACCUM_W_P(ACCUM_W), .PROD_W_P(PROD_W)
    ) u_dut (
        .clk(clk),             .rst_n(rst_n),
        .mm_wr_en(mm_wr_en),   .mm_wr_addr(mm_wr_addr), .mm_wr_data(mm_wr_data),
        .active_n(active_n_in),.num_edges(num_edges_in),
        .wr_xy_idx(wr_xy_idx), .wr_xy_en(wr_xy_en),
        .wr_x(wr_x),           .wr_y(wr_y),
        .run(run),             .Nstep(Nstep_in),
        .a0_fp(a0_fp_in),      .dt_fp(dt_fp_in), .c0_fp(c0_fp_in),
        .x_out(x_out),         .y_out(y_out),
        .signs_out(signs_out), .step_done(step_done),
        .schedule_done(schedule_done), .scan_done(scan_done)
    );

    // ── Edge list (0-indexed, matching G_simple.txt) ──────────────────────────
    // Format: {row, col, weight=1}
    typedef struct { int r, c, w; } edge_t;
    edge_t EDGES [0:NUM_EDGES-1] = '{
        '{0,1,-1}, '{0,2,-1}, '{1,3,-1}, '{2,3,-1},
        '{3,4,-1}, '{3,5,-1}, '{4,6,-1}, '{5,6,-1},
        '{6,7,-1}, '{0,7,-1}
    };

    // ── Initial x values (small, alternating sign to break symmetry) ──────────
    // Odd nodes start slightly negative, even nodes slightly positive
    // so the algorithm has a small hint toward the correct partition
    logic signed [XY_W-1:0] X_INIT [0:NUM_NODES-1] = '{
        16'sd50,    // node 0: slight +  (partition A)
       -16'sd30,    // node 1: slight -  (partition B)
       -16'sd20,    // node 2: slight -  (partition B)
        16'sd40,    // node 3: slight +  (partition A)
       -16'sd35,    // node 4: slight -  (partition B)
       -16'sd25,    // node 5: slight -  (partition B)
        16'sd45,    // node 6: slight +  (partition A)
       -16'sd15     // node 7: slight -  (partition B) -- contested, watch this one
    };

    // ── Tasks ─────────────────────────────────────────────────────────────────
    task automatic mm_write(
        input logic [MAIN_ADDR_W-1:0] addr,
        input logic [MAIN_WORD_W-1:0] data
    );
        @(posedge clk);
        mm_wr_en = 1'b1; mm_wr_addr = addr; mm_wr_data = data;
        @(posedge clk);
        mm_wr_en = 1'b0;
    endtask

    function automatic logic [MAIN_WORD_W-1:0] pack_edge(
        input int r, c, w
    );
        return MAIN_WORD_W'({COL_IDX_W'(r), COL_IDX_W'(c), IC_BITS'(signed'(w))});
    endfunction

    task automatic wait_high(ref logic sig, input int timeout, input string lbl);
        int cnt = 0;
        while (!sig) begin
            @(posedge clk);
            if (++cnt > timeout) $fatal(1, "TIMEOUT: %s", lbl);
        end
    endtask

    // ── Stimulus ──────────────────────────────────────────────────────────────
    initial begin
        // Full VCD dump -- this is a tiny run so the file will be small
        $dumpfile("tb_simple.vcd");
        $dumpvars(0, tb_simple);

        $display("");
        $display("====================================================");
        $display(" Simple 8-node MaxCut waveform test");
        $display(" Nodes=8  Edges=10  Nstep=%0d", NSTEP);
        $display(" Optimal cut = 10 (all edges, bipartite graph)");
        $display(" Optimal partition: A={0,3,6}  B={1,2,4,5,7}");
        $display("====================================================");

        // ── Reset ─────────────────────────────────────────────────────────────
        repeat(4) @(posedge clk);
        rst_n = 1'b1;
        repeat(2) @(posedge clk);
        $display("[%0t] Reset released", $time);

        // ── Phase 1: Load edges ───────────────────────────────────────────────
        $display("[%0t] Phase 1: Loading %0d edges", $time, NUM_EDGES);
        for (int e = 0; e < NUM_EDGES; e++) begin
            mm_write(MAIN_ADDR_W'(MAIN_EDGE_BASE + e),
                     pack_edge(EDGES[e].r, EDGES[e].c, EDGES[e].w));
            $display("  [%0t] edge[%0d]: %0d--%0d w=%0d",
                     $time, e, EDGES[e].r, EDGES[e].c, EDGES[e].w);
        end

        // ── Phase 2: Init x/y ─────────────────────────────────────────────────
        $display("[%0t] Phase 2: Initialising x", $time);
        for (int i = 0; i < NUM_NODES; i++) begin
            @(posedge clk);
            wr_xy_en  = 1'b1;
            wr_xy_idx = $clog2(N)'(i);
            wr_x      = X_INIT[i];
            wr_y      = '0;
            @(posedge clk);
            wr_xy_en  = 1'b0;
            $display("  [%0t] x[%0d] = %0d", $time, i, $signed(X_INIT[i]));
        end

        // ── Phase 3: Configure and run ────────────────────────────────────────
        $display("[%0t] Phase 3: Setting active_n=%0d num_edges=%0d",
                 $time, NUM_NODES, NUM_EDGES);
        active_n_in  = $clog2(N)'(NUM_NODES);
        num_edges_in = CSR_ADDR_W'(NUM_EDGES);

        @(posedge clk);
        $display("[%0t] Asserting run=1 (scan will start automatically)", $time);
        run = 1'b1;

        // ── Wait for scan ──────────────────────────────────────────────────────
        wait_high(scan_done, 10000, "scan_done");
        $display("[%0t] scan_done -- CSR store built", $time);

        // ── Wait for solver ────────────────────────────────────────────────────
        $display("[%0t] Running %0d dSB steps...", $time, NSTEP);
        wait_high(schedule_done, 1_000_000, "schedule_done");
        @(posedge clk);
        run = 1'b0;
        $display("[%0t] schedule_done", $time);

        // ── Phase 4: Print and check result ───────────────────────────────────
        $display("");
        $display("====================================================");
        $display(" FINAL STATE");
        $display("====================================================");
        $display(" signs_out = %08b", signs_out);
        $display(" (bit 0 = node 0, bit 7 = node 7)");
        $display(" Expected optimal: signs = 0b01001001 = 8'b01001001");
        $display("   nodes {0,3,6} = 1 (+), nodes {1,2,4,5,7} = 0 (-)");
        $display("");

        for (int i = 0; i < NUM_NODES; i++) begin
            automatic logic signed [XY_W-1:0] xi;
            xi = x_out[i*XY_W +: XY_W];
            $display("  node %0d: x=%6d  sign=%0b  partition=%s",
                     i, $signed(xi), signs_out[i],
                     signs_out[i] ? "A (+)" : "B (-)");
        end

        // Compute cut
        begin
            automatic int cut = 0;
            for (int e = 0; e < NUM_EDGES; e++) begin
                if (signs_out[EDGES[e].r] !== signs_out[EDGES[e].c]) begin
                    cut++;
                    $display("  edge %0d--%0d: CUT", EDGES[e].r, EDGES[e].c);
                end else begin
                    $display("  edge %0d--%0d: not cut (both in partition %s)",
                             EDGES[e].r, EDGES[e].c,
                             signs_out[EDGES[e].r] ? "A" : "B");
                end
            end
            $display("");
            $display("  Achieved cut = %0d / 10 (optimal = 10)", cut);
            if (cut == 10)
                $display("  PASS -- optimal solution found!");
            else if (cut >= 8)
                $display("  GOOD -- near-optimal (%0d/10)", cut);
            else
                $display("  BELOW EXPECTATION -- cut=%0d, check waveform", cut);
        end

        $display("====================================================");
        $display("");
        $display("Waveform saved to tb_simple.vcd");
        $display("Open in SimVision and probe these signals first:");
        $display("  tb_simple.clk");
        $display("  tb_simple.rst_n");
        $display("  tb_simple.run");
        $display("  tb_simple.scan_done");
        $display("  tb_simple.schedule_done");
        $display("  tb_simple.step_done");
        $display("  tb_simple.signs_out");
        $display("  tb_simple.u_dut.u_scan.*    (scanner FSM state)");
        $display("  tb_simple.u_dut.u_core.*    (SLC state, Jx_i value)");
        $display("  tb_simple.u_dut.x_reg[0..7] (oscillator positions)");

        repeat(4) @(posedge clk);
        $finish;
    end

    // ── Watchdog ──────────────────────────────────────────────────────────────
    initial begin
        #5_000_000;
        $fatal(1, "GLOBAL TIMEOUT");
    end

    // ── Per-step monitor: print sign vector and step number after each step ───
    int step_count = 0;
    always @(posedge clk) begin
        if (step_done) begin
            step_count++;
            $display("[%0t] step %0d done  signs=%08b  (A={0,3,6} would be 01001001)",
                     $time, step_count, signs_out);
        end
    end

endmodule

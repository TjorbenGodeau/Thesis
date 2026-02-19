// filepath: c:/Users/tjorb/Documents/Thesis/rtl/testbench.sv
// Testbench for GbSB module

`timescale 1ns / 1ps

module testbench;

    // Parameters
    parameter int N = 8;
    parameter int DATA_WIDTH = 32;
    parameter int FRAC_WIDTH = 16;
    parameter bit PER_SPIN = 1'b0;
    parameter real CLK_PERIOD = 10.0; // 10ns clock

    // Signals
    logic clk;
    logic rst_n;
    logic start;
    logic [31:0] M;
    logic done;
    logic [31:0] step_count;
    
    logic init_valid;
    logic [N-1:0][DATA_WIDTH-1:0] x0;
    logic [N-1:0][DATA_WIDTH-1:0] y0;
    logic [N-1:0][DATA_WIDTH-1:0] p0;
    
    logic [N-1:0][N-1:0][DATA_WIDTH-1:0] J;
    
    logic [DATA_WIDTH-1:0] dt;
    logic [DATA_WIDTH-1:0] A;
    
    logic valid_out;
    logic [N-1:0][DATA_WIDTH-1:0] x;
    logic [N-1:0][DATA_WIDTH-1:0] y;
    logic [DATA_WIDTH-1:0] p_global;
    logic [N-1:0][DATA_WIDTH-1:0] p_individual;
    
    logic [N-1:0] spins;
    logic [DATA_WIDTH-1:0] energy;
    logic [DATA_WIDTH-1:0] hamiltonian;
    
    logic hit_boundary;
    logic [N-1:0] boundary_mask;

    // Instantiate DUT
    GbSB #(
        .N(N),
        .DATA_WIDTH(DATA_WIDTH),
        .FRAC_WIDTH(FRAC_WIDTH),
        .PER_SPIN(PER_SPIN)
    ) dut (
        .clk(clk),
        .rst_n(rst_n),
        .start(start),
        .M(M),
        .done(done),
        .step_count(step_count),
        .init_valid(init_valid),
        .x0(x0),
        .y0(y0),
        .p0(p0),
        .J(J),
        .dt(dt),
        .A(A),
        .valid_out(valid_out),
        .x(x),
        .y(y),
        .p_global(p_global),
        .p_individual(p_individual),
        .spins(spins),
        .energy(energy),
        .hamiltonian(hamiltonian),
        .hit_boundary(hit_boundary),
        .boundary_mask(boundary_mask)
    );

    // Clock generation
    initial begin
        clk = 1'b0;
        forever #(CLK_PERIOD/2.0) clk = ~clk;
    end

    // Test stimulus
    initial begin
        // Initialize signals
        rst_n = 1'b0;
        start = 1'b0;
        M = 32'd50;
        init_valid = 1'b0;
        dt = 32'h0001_0000; // 1.0 in fixed-point (FRAC_WIDTH=16)
        A = 32'h0000_0000;  // 0.0 in fixed-point
        
        // Reset
        #(CLK_PERIOD * 5) rst_n = 1'b1;
        
        // Initialize coupling matrix J (simple diagonal)
        for (int i = 0; i < N; i++) begin
            for (int j = 0; j < N; j++) begin
                if (i == j) begin
                    J[i][j] = 32'h0000_0000;
                end else begin
                    J[i][j] = 32'h0000_8000; // 0.5 in fixed-point
                end
            end
        end
        
        // Initialize positions x0 in range [-1.0, 1.0]
        for (int i = 0; i < N; i++) begin
            x0[i] = 32'h0000_0000; // Start at 0
            y0[i] = 32'h0000_0000; // Zero momentum
            p0[i] = 32'h0001_0000; // p = 1.0
        end
        
        // Set initial conditions
        init_valid = 1'b1;
        #(CLK_PERIOD);
        init_valid = 1'b0;
        
        // Start computation
        start = 1'b1;
        #(CLK_PERIOD);
        start = 1'b0;
        
        // Wait for completion
        wait(done);
        #(CLK_PERIOD * 5);
        
        // Display results
        $display("=== Simulation Complete ===");
        $display("Final Step Count: %0d", step_count);
        $display("Final Energy: 0x%08x", energy);
        $display("Final Hamiltonian: 0x%08x", hamiltonian);
        $display("Spins: %b", spins);
        $display("Boundary Hit: %b", hit_boundary);
        
        #(CLK_PERIOD * 10);
        $finish;
    end

    // Monitor outputs
    initial begin
        $monitor("Time=%0t | step=%0d | valid=%b | x[0]=0x%08x | y[0]=0x%08x | boundary=%b",
                 $time, step_count, valid_out, x[0], y[0], hit_boundary);
    end

endmodule
// tb_schedule.sv
`timescale 1ns/1ps

module tb_schedule;

    import dsb_pkg::*;

    // Parameters for test
    localparam A_BITS_T  = 16;
    localparam STEP_W_T  = 12;

    // Inputs
    logic [STEP_W_T-1:0] m;
    logic [STEP_W_T-1:0] Nstep;
    logic [A_BITS_T-1:0] a0_fp;
    // Output
    logic [A_BITS_T-1:0] a_m;

    // DUT
    schedule #(
        .A_BITS_P(A_BITS_T),
        .STEP_W_P(STEP_W_T)
    ) u_dut (
        .m(m),
        .Nstep(Nstep),
        .a0_fp(a0_fp),
        .a_m(a_m)
    );

    // Test control
    int error_count = 0;

    // Helper to compute expected value using integer arithmetic (no rounding)
    function automatic logic [A_BITS_T-1:0] expected_a_m(
        input logic [STEP_W_T-1:0] m,
        input logic [STEP_W_T-1:0] Nstep,
        input logic [A_BITS_T-1:0] a0_fp
    );
        logic [A_BITS_T+STEP_W_T-1:0] numerator;
        logic [A_BITS_T-1:0] result;
        if (m >= Nstep) begin
            result = '0;
        end else begin
            numerator = a0_fp * (Nstep - m);
            result = numerator / Nstep;  // integer division truncates
        end
        return result;
    endfunction

    // Task to check a test case
    task check_case(input string label,
                    input logic [STEP_W_T-1:0] m_in,
                    input logic [STEP_W_T-1:0] Nstep_in,
                    input logic [A_BITS_T-1:0] a0_in);
        logic [A_BITS_T-1:0] expected;
        expected = expected_a_m(m_in, Nstep_in, a0_in);
        // Apply inputs (combinational, so no clock needed)
        m = m_in; Nstep = Nstep_in; a0_fp = a0_in;
        #1; // let combinational settle
        if (a_m !== expected) begin
            $error("%s: m=%d, Nstep=%d, a0=%d: expected %d, got %d",
                   label, m_in, Nstep_in, a0_in, expected, a_m);
            error_count++;
        end else begin
            $display("%s: OK (a_m=%d)", label, a_m);
        end
    endtask

    // Test sequence
    initial begin
        $display("=== Starting schedule module test ===");

        // Test case 1: Nstep=20, a0=16384 (1.0 in Q0.14)
        check_case("Test 1a: m=0",  0, 20, 16384);
        check_case("Test 1b: m=5",  5, 20, 16384);
        check_case("Test 1c: m=10", 10, 20, 16384);
        check_case("Test 1d: m=15", 15, 20, 16384);
        check_case("Test 1e: m=19", 19, 20, 16384);
        check_case("Test 1f: m=20 (>=Nstep)", 20, 20, 16384);

        // Test case 2: Nstep=100, a0=8192 (0.5)
        check_case("Test 2a: m=0",  0, 100, 8192);
        check_case("Test 2b: m=50", 50, 100, 8192);
        check_case("Test 2c: m=99", 99, 100, 8192);
        check_case("Test 2d: m=100", 100, 100, 8192);

        // Test case 3: a0=0
        check_case("Test 3a: m=0", 0, 20, 0);
        check_case("Test 3b: m=10", 10, 20, 0);

        // Test case 4: Nstep=1 (only one step)
        check_case("Test 4a: m=0, Nstep=1", 0, 1, 16384);
        check_case("Test 4b: m=1, Nstep=1", 1, 1, 16384); // should be 0

        // Test case 5: Nstep large, a0 max
        check_case("Test 5a: m=0, Nstep=4095, a0=65535", 0, 4095, 65535);
        check_case("Test 5b: m=2047, Nstep=4095, a0=65535", 2047, 4095, 65535);

        // Test case 6: Check rounding/truncation - integer division truncates toward zero
        // a0=1000, Nstep=10, m=3: numerator=7000, /10=700 exactly
        check_case("Test 6a: exact division", 3, 10, 1000);
        // a0=999, Nstep=10, m=3: numerator=6993, /10=699 (truncated)
        check_case("Test 6b: truncation", 3, 10, 999);

        // Summary
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
        $dumpfile("tb_schedule.vcd");
        $dumpvars(0, tb_schedule);
    end

endmodule
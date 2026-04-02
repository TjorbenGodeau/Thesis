module xnor_phase1
    import dsb_pkg::*;
    #(
        parameter int unsigned N_P = N,
        parameter int unsigned IC_BITS_P = IC_BITS
    )(
        input logic clk,
        input logic capture,                        // strobe: latch on posedge
        dsb_tile_if.ph1_in tile_bus,                
        output logic [N_P*IC_BITS_P-1:0] xnor_J,
        output logic [N_P-1:0] sign_eq              // 1 = sale sign as x_i
    );

    always_ff @(posedge clk) begin
        if (capture) begin
            xnor_J <= tile_bus.rbl_J;
            sign_eq <= tile_bus.rbl_signs;
        end
    end
endmodule
module compute_tile
    import dsb_pkg::*;
#(
    parameter int unsigned N_P = N,
    parameter int unsigned IC_BITS_P = IC_BITS
)(
    input logic clk,
    input logic precharge,
    // Normal write port
    input logic wwl,
    input logic [N_P*IC_BITS_P-1:0] wbl_J,
    input logic [N_P-1:0] wbl_signs,
    // Compute port: target spin drives RWL
    input logic sign_xi,
    // Outputs via interface
    dsb_tile_if.tile_out tile_bus
);

    // J bitcells: N_P oscillators x IC_BITS_P bits each
    for (genvar n = 0; n < N_P; n++) begin : gen_osc
        for (genvar b = 0; b < IC_BITS_P; b++) begin : gen_bit
        dsb_8t_bitcell u_J (
            .clk (clk),
            .wwl (wwl),
            .wbl (wbl_J[n*IC_BITS_P + b]),
            .rwl (sign_xi),
            .precharge (precharge),
            .rbl_out (tile_bus.rbl_J[n*IC_BITS_P + b])
        );
        end
    end

    // Sign bitcells: one per oscillator
    for (genvar n = 0; n < N_P; n++) begin : gen_sign
        dsb_8t_bitcell u_s (
            .clk (clk),
            .wwl (wwl),
            .wbl (wbl_signs[n]),
            .rwl (sign_xi),
            .precharge (precharge),
            .rbl_out (tile_bus.rbl_signs[n])
        );
    end
endmodule
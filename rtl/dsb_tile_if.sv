interface dsb_tile_if
    import dsb_pkg::*;
    #(
        parameter int unsigned N_P = N,
        parameter int unsigned IC_BITS_P = IC_BITS
    );

    logic [N_P*IC_BITS_P-1:0] rbl_J;
    logic [N_P-1:0] rbl_signs;

    // Tile drives these signals toward the pashe-1 capture register
    modport tile_out (output rbl_J, output rbl_signs);
    modport ph1_in (input rbl_J, input rbl_signs);

endinterface
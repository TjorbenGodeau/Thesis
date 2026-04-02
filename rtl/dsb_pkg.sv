package dsb_pkg;
    // --- Problem size ---
    parameter int unsigned N = 8;       // number of dSB oscillators
    parameter int unsigned IC_BITS = 7; // J coefficient bit width

    // --- Fixed-point format for x and y ---
    parameter int unsigned XY_INT = 2;
    parameter int unsigned XY_FRAC = 14;
    parameter int unsigned XY_W = XY_INT + XY_FRAC;

    // --- Schedule coefficient a(m) ---
    parameter int unsigned A_BITS = 16;

    // --- Step counter ---
    parameter int unsigned STEP_W = 12;     // supports Nstep up to 2^STEP_W - 1

    // --- JX accumulator (must hold N * max|J|) ---
    parameter int unsigned ACCUM_W = 24;

    // --- a*x product width ---
    parameter int unsigned PROD_W = XY_W + A_BITS;

    // --- Fixed-point unity for x ---
    parameter logic signed [XY_W-1:0] ONE_FP = XY_W'(1 << XY_FRAC);

    // --- Core FSM states ---
    typedef enum logic [2:0] {
        S_IDLE = 3'd0,
        S_PRECHARGE = 3'd1,
        S_PHASE1 = 3'd2,
        S_PHASE2 = 3'd3,
        S_PHASE3 = 3'd4
    } core_state_t;

    // --- Array controller FSM states ---
    typedef enum logic [1:0] {
        CTRL_IDLE = 2'd0,
        CTRL_LOAD = 2'd1,
        CTRL_RUN = 2'd2
    } ctrl_state_t;
endpackage
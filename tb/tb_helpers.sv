package tb_helpers;
  import dsb_pkg::*;
 
  // Convert a real value to Q(XY_INT).(XY_FRAC) fixed-point
  function automatic logic signed [XY_W-1:0] to_fp (input real v);
    return XY_W'($rtoi(v * real'(1 << XY_FRAC)));
  endfunction
 
  // Convert fixed-point back to real
  function automatic real from_fp (input logic signed [XY_W-1:0] v);
    return real'($signed(v)) / real'(1 << XY_FRAC);
  endfunction
 
  // Sign of a fixed-point x value: 1 = positive/zero, 0 = negative
  function automatic logic sign_of (input logic signed [XY_W-1:0] v);
    return ~v[XY_W-1];
  endfunction
 
  // Encode a real J value to IC_BITS 2's complement
  // J range assumed [-1,+1], scaled so +1.0 → 8'h01 (integer encoding)
  function automatic logic [IC_BITS-1:0] J_encode (input int j_int);
    return IC_BITS'(j_int);
  endfunction
 
  // Tolerance for fixed-point comparison (1 LSB)
  parameter real TOL = 1.0 / real'(1 << XY_FRAC);
 
endpackage : tb_helpers
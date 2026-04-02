// =============================================================================
// TB 2 — dsb_compute_tile
// Tests:
//   T2.1  Write J row all-ones, signs all-one; sign_xi=1 → all XNOR J bits=1,
//         all sign_eq=1
//   T2.2  sign_xi=0 → XNOR with stored 1s produces 0, sign_eq becomes 0
//   T2.3  Mixed J: alternating 0/1 bits, verify per-bit XNOR
//   T2.4  Precharge releases correctly after one cycle
// =============================================================================
module tb_compute_tile;
  import dsb_pkg::*;
 
  localparam int unsigned N_T      = 4;
  localparam int unsigned IC_T     = 4;
 
  logic clk = 0;
  always #5 clk = ~clk;
 
  logic                   precharge, wwl, sign_xi;
  logic [N_T*IC_T-1:0]   wbl_J;
  logic [N_T-1:0]         wbl_signs;
 
  dsb_tile_if #(.N_P(N_T), .IC_BITS_P(IC_T)) tile_bus ();
 
  dsb_compute_tile #(.N_P(N_T), .IC_BITS_P(IC_T)) dut (
    .clk       (clk),
    .precharge (precharge),
    .wwl       (wwl),
    .wbl_J     (wbl_J),
    .wbl_signs (wbl_signs),
    .sign_xi   (sign_xi),
    .tile_bus  (tile_bus.tile_out)
  );
 
  int pass_cnt = 0, fail_cnt = 0;
 
  task automatic check_bits(
    input string             label,
    input logic [N_T*IC_T-1:0] got_J,
    input logic [N_T*IC_T-1:0] exp_J,
    input logic [N_T-1:0]      got_s,
    input logic [N_T-1:0]      exp_s
  );
    if (got_J === exp_J && got_s === exp_s) begin
      $display("  PASS [%s]", label);
      pass_cnt++;
    end else begin
      $error("  FAIL [%s] J: got=%0b exp=%0b  signs: got=%0b exp=%0b",
             label, got_J, exp_J, got_s, exp_s);
      fail_cnt++;
    end
  endtask
 
  task automatic write_tile(
    input logic [N_T*IC_T-1:0] j_data,
    input logic [N_T-1:0]       s_data
  );
    @(posedge clk); #1;
    wwl = 1; wbl_J = j_data; wbl_signs = s_data; precharge = 0; sign_xi = 0;
    @(posedge clk); #1;
    wwl = 0;
  endtask
 
  initial begin
    $dumpfile("tb_compute_tile.vcd");
    $dumpvars(0, tb_compute_tile);
    $display("=== TB2: dsb_compute_tile ===");
 
    wwl = 0; precharge = 0; sign_xi = 0; wbl_J = '0; wbl_signs = '0;
    repeat(2) @(posedge clk);
 
    // T2.1 All J bits=1, all signs=1; sign_xi=1 → XNOR(1,1)=1 everywhere
    write_tile({N_T*IC_T{1'b1}}, {N_T{1'b1}});
    @(posedge clk); #1; precharge = 1;
    @(posedge clk); #1; precharge = 0; sign_xi = 1;
    #2;
    check_bits("T2.1 allone_xi1",
               tile_bus.rbl_J,    {N_T*IC_T{1'b1}},
               tile_bus.rbl_signs, {N_T{1'b1}});
 
    // T2.2 All J=1, signs=1; sign_xi=0 → XNOR(1,0)=0 everywhere
    @(posedge clk); #1; sign_xi = 0;
    #2;
    check_bits("T2.2 allone_xi0",
               tile_bus.rbl_J,    {N_T*IC_T{1'b0}},
               tile_bus.rbl_signs, {N_T{1'b0}});
 
    // T2.3 Alternating J bits: 4'b1010 per oscillator; sign_xi=1
    // XNOR(1,1)=1, XNOR(0,1)=0 → per-nibble: 1010
    write_tile({N_T{4'b1010}}, {N_T{1'b1}});
    @(posedge clk); #1; precharge = 1;
    @(posedge clk); #1; precharge = 0; sign_xi = 1;
    #2;
    check_bits("T2.3 alt_bits_xi1",
               tile_bus.rbl_J,    {N_T{4'b1010}},
               tile_bus.rbl_signs, {N_T{1'b1}});
 
    // T2.4 Precharge releases → bits should follow XNOR again
    @(posedge clk); #1; precharge = 1;
    #1;
    // During precharge all rbl_J must be 1
    if (tile_bus.rbl_J === {N_T*IC_T{1'b1}})
      begin $display("  PASS [T2.4 precharge_high]"); pass_cnt++; end
    else
      begin $error("  FAIL [T2.4 precharge_high] rbl_J=%0b", tile_bus.rbl_J); fail_cnt++; end
    @(posedge clk); #1; precharge = 0; sign_xi = 1;
    #2;
    check_bits("T2.4 post_precharge",
               tile_bus.rbl_J,    {N_T{4'b1010}},
               tile_bus.rbl_signs, {N_T{1'b1}});
 
    $display("--- TB2 result: %0d PASS, %0d FAIL ---\n", pass_cnt, fail_cnt);
    #20 $finish;
  end
 
endmodule : tb_compute_tile
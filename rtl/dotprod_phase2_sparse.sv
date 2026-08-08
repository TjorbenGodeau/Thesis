// =============================================================================
// dotprod_phase2_sparse.sv
//
// Extended version of dotprod_phase2 that supports chunked (partial)
// accumulation for sparse rows that overflow one K_MAX compute tile.
//
// Key differences from the original dotprod_phase2
// ------------------------------------------------
//  1. N_P parameter is replaced by K_MAX_P (entries per chunk).
//     The tile only contains K_MAX_P valid entries per chunk, not N.
//
//  2. New input: accumulate
//       0 → fresh start: clear the accumulator before this chunk
//       1 → add to existing Jx_i (continuing from a previous chunk)
//
//  3. New input: valid_count  [K_CNT_W bits]
//       Number of valid entries in the current chunk (1..K_MAX_P).
//       Entries beyond valid_count are ignored (tile slots may be partially
//       filled on the last chunk of a row).
//
//  4. sign_xi is still required (for the global +1 correction).
//       On chunked rows it must be the SAME sign_xi for all chunks of row i
//       and the correction is applied ONLY on the LAST chunk
//       (caller sets last_chunk=1 on the final chunk).
//
//  5. Output Jx_i is now a registered accumulator that persists between
//     chunk calls (until accumulate=0 resets it).
//
// Interface unchanged otherwise: start/done handshake.
//
// Bit-width note
// --------------
//  With sparse entries the maximum |Jx_i| is still bounded by
//  N * max|J|  (same as the dense case).  ACCUM_W_P is unchanged.
//
// Timing
// ------
//  Same as original: done fires one cycle after start.
// =============================================================================

`timescale 1ns/1ps
module dotprod_phase2_sparse
    import dsb_pkg::*;
    import dsb_sparse_pkg::*;
#(
    parameter int unsigned K_MAX_P    = K_MAX,       // entries per chunk
    parameter int unsigned IC_BITS_P  = IC_BITS,
    parameter int unsigned ACCUM_W_P  = ACCUM_W,
    parameter int unsigned K_CNT_W    = $clog2(K_MAX + 1)
)(
    input  logic                         clk,

    // ── Handshake ─────────────────────────────────────────────────────────────
    input  logic                         start,       // pulse to trigger
    output logic                         done,        // pulses one cycle later

    // ── Chunk control ─────────────────────────────────────────────────────────
    input  logic                         accumulate,  // 0=fresh, 1=add to prev
    input  logic                         last_chunk,  // 1 on final chunk of row
    input  logic [K_CNT_W-1:0]           valid_count, // valid entries this chunk

    // ── Data from compute tile (same as original, but K_MAX_P wide) ───────────
    input  logic                         sign_xi,
    input  logic [K_MAX_P*IC_BITS_P-1:0] xnor_J,     // from xnor_phase1
    input  logic [K_MAX_P-1:0]           sign_eq,     // from xnor_phase1

    // ── Result ────────────────────────────────────────────────────────────────
    output logic signed [ACCUM_W_P-1:0]  Jx_i        // final or partial result
);

    // ── Persistent accumulator ────────────────────────────────────────────────
    // Holds the running sum across chunks for one oscillator row.
    logic signed [ACCUM_W_P-1:0] Jx_accum;

    always_ff @(posedge clk) begin
        done <= 1'b0;

        if (start) begin
            // ── Build the chunk partial sum ───────────────────────────────────
            automatic logic signed [ACCUM_W_P-1:0] chunk_sum;
            chunk_sum = '0;

            for (int k = 0; k < K_MAX_P; k++) begin
                // Only accumulate valid entries
                if (k < int'(valid_count)) begin
                    automatic logic [IC_BITS_P-1:0] xnor_word;
                    automatic logic [IC_BITS_P-1:0] corrected;
                    xnor_word = xnor_J[k*IC_BITS_P +: IC_BITS_P];
                    // sign-equality correction (same logic as original)
                    corrected  = sign_eq[k] ? xnor_word : ~xnor_word;
                    chunk_sum += ACCUM_W_P'(signed'(corrected));
                end
            end

            // ── Add to or replace accumulator, with sign correction ───────────
            // When sign(x_i) is negative the XNOR encoding adds an implicit
            // +1 per term (same correction as the original dotprod_phase2).
            // For chunked rows the correction is deferred to the last chunk,
            // where the caller supplies the TOTAL nnz of the row as valid_count
            // so we can correct for all accumulated terms at once.
            //
            // Base value: prev_accum + chunk_sum
            //   prev_accum = Jx_accum if accumulate, else 0
            begin
                automatic logic signed [ACCUM_W_P-1:0] base;
                base = (accumulate ? Jx_accum : '0) + chunk_sum;

                if (last_chunk && !sign_xi)
                    // valid_count == total nnz of the row (set by SLC)
                    Jx_accum <= base + ACCUM_W_P'(valid_count);
                else
                    Jx_accum <= base;
            end

            done <= 1'b1;
        end
    end

    // ── Drive output from accumulator ─────────────────────────────────────────
    assign Jx_i = Jx_accum;

endmodule

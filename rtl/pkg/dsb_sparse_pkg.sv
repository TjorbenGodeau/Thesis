// =============================================================================
// dsb_sparse_pkg.sv   (v3 — Gset-scale rewrite)
//
// Extension of dsb_pkg with parameters and types for the sparse memory
// subsystem.  Import BOTH packages where needed:
//   import dsb_pkg::*;
//   import dsb_sparse_pkg::*;
//
// ─────────────────────────────────────────────────────────────────────────────
// WHY THIS FILE CHANGED FROM THE PROTOTYPE VERSION
// ─────────────────────────────────────────────────────────────────────────────
// The original prototype addressed main memory densely: J[i][j] lived at
// address i*N+j, requiring N*N words.  That is fine for N=8 (64 words) but
// is fundamentally broken for Gset-scale problems:
//
//   N = 20,000  →  N*N = 400,000,000 words  →  800 MB at 16 bits/word
//
// This is absurd for a design whose entire point is exploiting sparsity —
// you would be allocating a dense matrix just to discover it is 99.97% zero.
//
// FIX: main_memory now stores the J matrix as a flat EDGE LIST:
//   edge[e] = {row_idx, col_idx, weight}   for e = 0 .. num_edges-1
// This is also a much more natural match for the actual Gset file format,
// which IS an edge list (see the format note below).
//
// The sparsity_scanner is renamed in spirit to an "edge list ingest" FSM:
// it no longer searches for non-zero entries (every edge in the list is by
// construction non-zero) — it just walks the edge list once and, for each
// edge (i,j,w), writes TWO CSR entries (i→j and j→i) because the J matrix
// is symmetric (undirected graph).
//
// ─────────────────────────────────────────────────────────────────────────────
// Gset file format (for reference — read by the testbench, not by hardware)
// ─────────────────────────────────────────────────────────────────────────────
//   Line 1:        <num_nodes> <num_edges>
//   Lines 2..end:  <node_i> <node_j> <weight>      (1-indexed nodes)
//
// ─────────────────────────────────────────────────────────────────────────────
// Memory map (word-addressable main memory) — v3
// ─────────────────────────────────────────────────────────────────────────────
//  [0 .. MAX_EDGES-1]                     : edge list, one edge per word
//                                            packed as {row, col, weight}
//  [MAX_EDGES .. MAX_EDGES+N-1]           : x_init[i]
//  [MAX_EDGES+N .. MAX_EDGES+2*N-1]       : y_init[i]
//
// ─────────────────────────────────────────────────────────────────────────────
// Sparse Compute Memory layout (per chunk in the 8T tile) — unchanged concept
// ─────────────────────────────────────────────────────────────────────────────
//  Each row holds K_MAX sparse entries.
//  Entry k  : col_index stored in plain SRAM  (COL_IDX_W bits)
//             J value   stored in 8T bitcells (IC_BITS bits)
//
// ─────────────────────────────────────────────────────────────────────────────
// CSR index store layout — unchanged concept, now sized for real Gset edges
// ─────────────────────────────────────────────────────────────────────────────
//  row_ptr[i]  : {start_addr, count}  for row i   (ROW_PTR_W wide)
//  entry[a]    : {col_index, edge_addr}            (CSR_ENTRY_W wide)
//                edge_addr points back into the edge-list region of
//                main_memory so the weight can be re-read when needed.
// =============================================================================

package dsb_sparse_pkg;

    import dsb_pkg::*;

    // -------------------------------------------------------------------------
    // Problem-size ceiling for Gset benchmarking
    // -------------------------------------------------------------------------
    // dsb_pkg::N must be set to the LARGEST instance you intend to run in a
    // given synthesis/simulation build (e.g. 20000 for G77/G81). Smaller
    // instances (e.g. 800-node G1..G67) simply use a subset 0..N_active-1;
    // the testbench drives a separate "active_n" so you do not need to
    // recompile dsb_pkg for every Gset file. See tb_gset_benchmark.sv.
    //
    // NOTE: dsb_pkg::N is a COMPILE-TIME parameter, so the widths below
    // (COL_IDX_W, etc.) are sized for whatever N you set in dsb_pkg at
    // build time.  Set N = 20000 to be able to run every Gset instance in
    // the same build (smaller graphs just leave the unused oscillators idle).
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Sparse compute memory geometry
    // -------------------------------------------------------------------------
    // Number of sparse-entry slots per physical compute-memory row/chunk.
    // Tune this to match your 8T array row width.  Larger K_MAX means fewer
    // chunks per row (less SLC overhead) at the cost of more 8T bitcells.
    parameter int unsigned K_MAX = 64;

    // -------------------------------------------------------------------------
    // Edge-list capacity
    // -------------------------------------------------------------------------
    // Maximum number of UNDIRECTED edges the design can hold.  Sized against
    // the largest Gset instance: G63/G64 (41,459 edges) and a safety margin.
    // Each undirected edge becomes 2 CSR entries (i->j and j->i).
    parameter int unsigned MAX_EDGES = 50_000;

    // Total directed CSR entries = 2 * MAX_EDGES (worst case, all edges used)
    parameter int unsigned MAX_NNZ   = 2 * MAX_EDGES;

    // -------------------------------------------------------------------------
    // Derived widths
    // -------------------------------------------------------------------------
    // Bits needed to index one of N oscillators (column/row index)
    parameter int unsigned COL_IDX_W = $clog2(N);

    // -------------------------------------------------------------------------
    // Main memory layout — EDGE LIST based (v3)
    // -------------------------------------------------------------------------
    // Each edge-list word packs {row_idx, col_idx, weight}.
    parameter int unsigned EDGE_ROW_W   = COL_IDX_W;
    parameter int unsigned EDGE_COL_W   = COL_IDX_W;
    parameter int unsigned EDGE_WT_W    = IC_BITS;
    parameter int unsigned EDGE_WORD_W  = EDGE_ROW_W + EDGE_COL_W + EDGE_WT_W;

    // Main memory must hold the WIDER of {edge word, XY_W word}, since x/y
    // init values are also stored there.
    parameter int unsigned MAIN_WORD_W  = (EDGE_WORD_W > XY_W) ? EDGE_WORD_W : XY_W;

    // Total word count: edge list + x_init + y_init
    parameter int unsigned MAIN_MEM_WORDS = MAX_EDGES + 2 * N;
    parameter int unsigned MAIN_ADDR_W    = $clog2(MAIN_MEM_WORDS);

    // Base addresses in main memory
    parameter int unsigned MAIN_EDGE_BASE = 0;
    parameter int unsigned MAIN_X_BASE    = MAX_EDGES;
    parameter int unsigned MAIN_Y_BASE    = MAX_EDGES + N;

    // -------------------------------------------------------------------------
    // CSR entry store address width
    // -------------------------------------------------------------------------
    parameter int unsigned CSR_ADDR_W  = $clog2(MAX_NNZ + 1);

    // One CSR entry: col_index + main_memory address of the edge word
    // (so the SLC can re-read the weight directly from the edge list)
    parameter int unsigned CSR_ENTRY_W = COL_IDX_W + MAIN_ADDR_W;

    // Row-pointer entry: start address into the CSR entry store + count.
    // Worst case one row holds ALL non-zeros (star graph), so size count
    // against MAX_NNZ for safety.
    parameter int unsigned ROW_COUNT_W = $clog2(MAX_NNZ + 1);
    parameter int unsigned ROW_PTR_W   = CSR_ADDR_W + ROW_COUNT_W;

    // -------------------------------------------------------------------------
    // Chunk counter: how many K_MAX-sized chunks fit in the worst-case row
    // -------------------------------------------------------------------------
    // Worst-case row length is bounded by MAX_NNZ (star-graph degree bound),
    // not by N, since Gset graphs are not guaranteed regular.
    parameter int unsigned MAX_ROW_NNZ = MAX_NNZ;
    parameter int unsigned MAX_CHUNKS  = (MAX_ROW_NNZ + K_MAX - 1) / K_MAX;
    parameter int unsigned CHUNK_CNT_W = $clog2(MAX_CHUNKS + 1);

    // -------------------------------------------------------------------------
    // FSM state types
    // -------------------------------------------------------------------------

    // Edge-list ingest FSM (was "sparsity_scanner" in the dense prototype;
    // kept the same module/signal names for interface compatibility, but the
    // states now reflect edge-list consumption instead of dense-grid scanning)
    typedef enum logic [2:0] {
        SCAN_IDLE     = 3'd0,
        SCAN_READ     = 3'd1,   // issue read of edge[e] from main memory
        SCAN_WAIT     = 3'd2,   // wait one cycle for BRAM output register
        SCAN_WR_FWD   = 3'd3,   // write CSR entry for i -> j
        SCAN_WR_REV   = 3'd4,   // write CSR entry for j -> i
        SCAN_NEXTEDGE = 3'd5,   // advance to next edge
        SCAN_FINALIZE = 3'd6,   // write all N row_ptr entries from row_count[]
        SCAN_DONE     = 3'd7
    } scan_state_t;

    // Sparse load controller FSM (unchanged from the prototype)
    typedef enum logic [3:0] {
        SLC_IDLE       = 4'd0,
        SLC_LOAD_PTR   = 4'd1,   // read row_ptr[i] from CSR store
        SLC_WAIT_PTR   = 4'd2,   // BRAM latency for row_ptr
        SLC_CHUNK_START= 4'd3,   // start filling one K_MAX chunk
        SLC_LOAD_ENTRY = 4'd4,   // read CSR entry[a]
        SLC_WAIT_ENTRY = 4'd5,   // BRAM latency for CSR entry
        SLC_READ_J     = 4'd6,   // read J value from main memory
        SLC_WAIT_J     = 4'd7,   // BRAM latency for J value
        SLC_WRITE_TILE = 4'd8,   // write (col, J) into compute memory slot
        SLC_NEXT_ENTRY = 4'd9,   // advance to next entry in this chunk
        SLC_TRIGGER_PH2= 4'd10,  // pulse dotprod start for this chunk
        SLC_WAIT_PH2   = 4'd11,  // wait for dotprod done
        SLC_NEXT_CHUNK = 4'd12,  // move to next chunk or finish
        SLC_DONE       = 4'd13
    } slc_state_t;

endpackage

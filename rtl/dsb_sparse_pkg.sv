// =============================================================================
// dsb_sparse_pkg.sv
//
// Extension of dsb_pkg with parameters and types for the sparse memory
// subsystem.  Import BOTH packages where needed:
//   import dsb_pkg::*;
//   import dsb_sparse_pkg::*;
//
// Memory map (word-addressable main memory)
// ------------------------------------------
//  [0 .. N*N - 1]              : J matrix  (J[i][j] at addr i*N + j)
//  [N*N .. N*N + N - 1]        : x_init    (x_init[i] at addr N*N + i)
//  [N*N + N .. N*N + 2*N - 1]  : y_init    (y_init[i] at addr N*N + N + i)
//
// Sparse Compute Memory layout (per row in the 8T tile)
// ------------------------------------------------------
//  Each row holds K_MAX sparse entries.
//  Entry k  : col_index stored in plain SRAM  (COL_IDX_W bits)
//             J value   stored in 8T bitcells (IC_BITS bits)
//
// CSR index store layout
// ----------------------
//  row_ptr[i]  : {start_addr, count}  for row i   (ROW_PTR_W wide)
//  entry[a]    : {col_index, main_mem_addr}        (CSR_ENTRY_W wide)
//                stored flat, row i's entries start at row_ptr[i].start
// =============================================================================

package dsb_sparse_pkg;

    import dsb_pkg::*;

    // -------------------------------------------------------------------------
    // Sparse compute memory geometry
    // -------------------------------------------------------------------------
    // Number of sparse-entry slots per physical compute-memory row.
    // Tune this to match your 8T array row width.
    parameter int unsigned K_MAX = 4;

    // Maximum non-zeros across the ENTIRE matrix (size the CSR entry store).
    // Default: assume at most N*(N-1)/2 unique edges for a dense upper triangle.
    parameter int unsigned MAX_NNZ = N * N;   // safe upper bound

    // -------------------------------------------------------------------------
    // Derived widths
    // -------------------------------------------------------------------------
    // Bits needed to index one of N oscillators (column index)
    parameter int unsigned COL_IDX_W = $clog2(N);

    // Bits needed to address the full main memory word space
    // Layout: J matrix (N*N words) + x_init (N words) + y_init (N words)
    parameter int unsigned MAIN_MEM_WORDS  = N * N + 2 * N;
    parameter int unsigned MAIN_ADDR_W     = $clog2(MAIN_MEM_WORDS);

    // Main memory word width: must hold the widest item.
    // J values are IC_BITS signed; x/y values are XY_W signed.
    // Store everything in XY_W-wide words (J values zero-extended on write).
    parameter int unsigned MAIN_WORD_W = XY_W;   // 16 bits default

    // CSR entry store address width
    parameter int unsigned CSR_ADDR_W  = $clog2(MAX_NNZ + 1);

    // One CSR entry: col_index + main_memory_address of that J coefficient
    parameter int unsigned CSR_ENTRY_W = COL_IDX_W + MAIN_ADDR_W;

    // Row-pointer entry: start address into the CSR entry store + count
    // count field is wide enough for MAX_NNZ entries in one row
    parameter int unsigned ROW_COUNT_W = $clog2(MAX_NNZ + 1);
    parameter int unsigned ROW_PTR_W   = CSR_ADDR_W + ROW_COUNT_W;

    // -------------------------------------------------------------------------
    // Chunk counter: how many K_MAX-sized chunks fit in the worst-case row
    // -------------------------------------------------------------------------
    parameter int unsigned MAX_CHUNKS  = (N + K_MAX - 1) / K_MAX;
    parameter int unsigned CHUNK_CNT_W = $clog2(MAX_CHUNKS + 1);

    // -------------------------------------------------------------------------
    // Base addresses in main memory
    // -------------------------------------------------------------------------
    parameter int unsigned MAIN_J_BASE = 0;
    parameter int unsigned MAIN_X_BASE = N * N;
    parameter int unsigned MAIN_Y_BASE = N * N + N;

    // -------------------------------------------------------------------------
    // FSM state types
    // -------------------------------------------------------------------------

    // Sparsity scanner FSM
    typedef enum logic [2:0] {
        SCAN_IDLE    = 3'd0,
        SCAN_READ    = 3'd1,   // issue read to main memory
        SCAN_WAIT    = 3'd2,   // wait one cycle for BRAM output register
        SCAN_EVAL    = 3'd3,   // check if non-zero, write CSR if so
        SCAN_NEXTCOL = 3'd4,   // advance column pointer
        SCAN_NEXTROW = 3'd5,   // close current row, write row_ptr, advance row
        SCAN_DONE    = 3'd6
    } scan_state_t;

    // Sparse load controller FSM
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

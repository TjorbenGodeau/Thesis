#!/bin/bash
# run_sparsity_scanner.sh

WORK_DIR="."   # set to your RTL directory if needed

xrun \
    -gui \
    -debug \
    -linedebug \
    $WORK_DIR/rtl/pkg/dsb_pkg.sv \
    $WORK_DIR/rtl/pkg/dsb_sparse_pkg.sv \
    $WORK_DIR/rtl/main_memory.sv \
    $WORK_DIR/rtl/csr_index_store.sv \
    $WORK_DIR/rtl/sparsity_scanner.sv \
    $WORK_DIR/tb/tb_sparsity_scanner.sv \
    -top tb_sparsity_scanner \
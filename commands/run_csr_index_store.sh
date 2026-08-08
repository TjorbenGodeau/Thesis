#!/bin/bash
# run_csr_index_store.sh

WORK_DIR="."   # set to your RTL directory if needed

xrun \
    -gui \
    -debug \
    -linedebug \
    $WORK_DIR/rtl/pkg/dsb_pkg.sv \
    $WORK_DIR/rtl/pkg/dsb_sparse_pkg.sv \
    $WORK_DIR/rtl/csr_index_store.sv \
    $WORK_DIR/tb/tb_csr_index_store.sv \
    -top tb_csr_index_store \
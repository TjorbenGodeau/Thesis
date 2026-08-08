#!/bin/bash
# run_sparse_load_controller.sh

WORK_DIR="."   # set to your RTL directory if needed

xrun \
    -gui \
    -debug \
    -linedebug \
    $WORK_DIR/rtl/pkg/dsb_pkg.sv \
    $WORK_DIR/rtl/pkg/dsb_sparse_pkg.sv \
    $WORK_DIR/rtl/main_memory.sv \
    $WORK_DIR/rtl/csr_index_store.sv \
    $WORK_DIR/rtl/sparse_compute_memory.sv \
    $WORK_DIR/rtl/8t_bitcell.sv \
    $WORK_DIR/rtl/dotprod_phase2_sparse.sv \
    $WORK_DIR/rtl/sparse_load_controller.sv \
    $WORK_DIR/tb/tb_sparse_load_controller.sv \
    -top tb_sparse_load_controller \
#!/bin/bash
# run_sparse_compute_memory.sh

# Set WORK_DIR to your RTL directory (where 8t_bitcell.sv and sparse_compute_memory.sv reside)
WORK_DIR="."

xrun \
    -gui \
    -debug \
    -linedebug \
    $WORK_DIR/rtl/pkg/dsb_pkg.sv \
    $WORK_DIR/rtl/pkg/dsb_sparse_pkg.sv \
    $WORK_DIR/rtl/8t_bitcell.sv \
    $WORK_DIR/rtl/sparse_compute_memory.sv \
    $WORK_DIR/tb/tb_sparse_compute_memory.sv \
    -top tb_sparse_compute_memory \
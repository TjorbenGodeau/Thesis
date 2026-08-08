#!/bin/bash
# run_dotprod_phase2_sparse.sh

WORK_DIR="."   # set to your RTL directory if needed

xrun \
    -gui \
    -debug \
    -linedebug \
    $WORK_DIR/rtl/pkg/dsb_pkg.sv \
    $WORK_DIR/rtl/pkg/dsb_sparse_pkg.sv \
    $WORK_DIR/rtl/dotprod_phase2_sparse.sv \
    $WORK_DIR/tb/tb_dotprod_phase2_sparse.sv \
    -top tb_dotprod_phase2_sparse \
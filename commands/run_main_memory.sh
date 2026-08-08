#!/bin/bash
# run_main_memory.sh

WORK_DIR="."   # set to your RTL directory if needed

xrun \
    -gui \
    -debug \
    -linedebug \
    $WORK_DIR/rtl/pkg/dsb_pkg.sv \
    $WORK_DIR/rtl/pkg/dsb_sparse_pkg.sv \
    $WORK_DIR/rtl/main_memory.sv \
    $WORK_DIR/tb/tb_main_memory.sv \
    -top tb_main_memory \
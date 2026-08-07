#!/bin/bash
# run_8t_bitcell_array.sh

WORK_DIR="."   # set to your RTL directory

xrun \
    -gui \
    -debug \
    -linedebug \
    $WORK_DIR/rtl/pkg/dsb_pkg.sv \
    $WORK_DIR/rtl/pkg/dsb_sparse_pkg.sv \
    $WORK_DIR/rtl/8t_bitcell.sv \
    $WORK_DIR/tb/tb_8t_bitcell_array.sv \
    -top tb_8t_bitcell_array \
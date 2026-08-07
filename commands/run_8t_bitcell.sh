#!/bin/bash
# run_8t_bitcell.sh

WORK_DIR="."   # set to your RTL directory if needed

xrun \
    -gui \
    -debug \
    -linedebug \
    $WORK_DIR/rtl/8t_bitcell.sv \
    $WORK_DIR/tb/tb_8t_bitcell.sv \
    -top tb_8t_bitcell \
    -exit
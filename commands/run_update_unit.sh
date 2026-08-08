#!/bin/bash
# run_update_unit.sh

WORK_DIR="."   # set to your RTL directory if needed

xrun \
    -gui \
    -debug \
    -linedebug \
    $WORK_DIR/rtl/pkg/dsb_pkg.sv \
    $WORK_DIR/rtl/update_unit.sv \
    $WORK_DIR/tb/tb_update_unit.sv \
    -top tb_update_unit \
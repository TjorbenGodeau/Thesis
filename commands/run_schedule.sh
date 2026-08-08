#!/bin/bash
# run_schedule.sh

WORK_DIR="."

xrun \
    -gui \
    -debug \
    -linedebug \
    $WORK_DIR/rtl/pkg/dsb_pkg.sv \
    $WORK_DIR/rtl/schedule.sv \
    $WORK_DIR/tb/tb_schedule.sv \
    -top tb_schedule \
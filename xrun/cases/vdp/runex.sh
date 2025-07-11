#!/usr/bin/env bash

NPROC=1
OUTPUT_FILE="out.txt"
FEXEC="vdp"

if [ ! -d "./datafiles" ]; then
    mkdir -p "./datafiles"
fi


# {PROJECT_DIR} will be replaced by the python script
BUILD_DIR="/home/khara/projects/petsc-ts/build"
cd $BUILD_DIR || exit; make ${FEXEC}; cd - || exit


mpirun -n ${NPROC} ${BUILD_DIR}/${FEXEC} > ${OUTPUT_FILE} 2>&1 \
-options_view \
-ts_monitor \
-ts_dt 0.001 \
-ts_max_time 200 \
-ts_save_trajectory
#-ts_type bdf \
#-ts_bdf_order 2 \
#-ts_adapt_type none \


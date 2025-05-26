#!/usr/bin/env bash

NPROC=1
OUTPUT_FILE="out.txt"
FEXEC="vanderpol_bvp"

if [ ! -d "./datafiles" ]; then
    mkdir -p "./datafiles"
fi

# {PROJECT_DIR} will be replaced by the python script
BUILD_DIR="/home/khara/projects/petsc-ts/build"
cd $BUILD_DIR || exit; make; cd - || exit

mpirun -n ${NPROC} ${BUILD_DIR}/${FEXEC} > ${OUTPUT_FILE} 2>&1 \
-options_view \
-da_refine 4 \
-snes_monitor \
-snes_rtol 1e-10 \
-snes_max_it 500 \
-snes_converged_reason \
#-ksp_monitor \
#-mat_view :matrix.m:ascii_matlab
# -da_refine_x 2 \
# -da_refine_y 2


python plot_sol.py
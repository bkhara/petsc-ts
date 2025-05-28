#!/usr/bin/env bash

NPROC=1
OUTPUT_FILE="out.txt"
FEXEC="vdp_periodic_ptc"

if [ ! -d "./datafiles" ]; then
    mkdir -p "./datafiles"
fi


# {PROJECT_DIR} will be replaced by the python script
BUILD_DIR="/home/khara/projects/petsc-ts/build"
cd $BUILD_DIR || exit; make ${FEXEC}; cd - || exit

mpirun -n ${NPROC} ${BUILD_DIR}/${FEXEC} > ${OUTPUT_FILE} 2>&1 \
-options_view \
-da_refine 8 \
-ts_type pseudo \
-ts_monitor_pseudo \
-ts_pseudo_frtol 1e-12 \
-ts_dt 1e-3 \
# -ts_monitor_solution_vtk datafiles/blah-%03D.vts \
# -ts_pseudo_increment 2 \
# -ts_pseudo_increment_dt_from_initial_dt
# -ts_max_time 200 \
# -ts_monitor \
#-ts_save_trajectory
#-ts_bdf_order 2 \
#-ts_adapt_type none \

python plot_sol.py
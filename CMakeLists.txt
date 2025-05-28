cmake_minimum_required(VERSION 3.20)
project(petsc_ts)

set(CMAKE_CXX_STANDARD 20)
# Add the cmake_modules directory to the CMake module search path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules")

if (NOT TARGET talyfem)
    find_package(talyfem REQUIRED)
endif()
find_package(P4EST REQUIRED)

add_executable(nptel_ode_ts src/nptel_examples/ode.cpp)
add_executable(nptel_poisson src/nptel_examples/poisson.cpp)
add_executable(nptel_heat2d_ts src/nptel_examples/heat2d.cpp)
add_executable(nptel_turing_patterns src/nptel_examples/turing_patterns.cpp)
add_executable(petsc_ex45 src/petsc_examples/ex45.c)
add_executable(p4est_simple src/p4est_simple.cpp)

# add the ./include directory to the include search path so we can #include
# the header files we put in include/
target_include_directories(nptel_ode_ts PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(nptel_poisson PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(nptel_heat2d_ts PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(nptel_turing_patterns PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(petsc_ex45 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(p4est_simple PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# this will also copy all the compiler definitions, etc. that TalyFEM was
# compiled with (that were marked as PUBLIC) - e.g. this will automatically
# add "-DENABLE_HDF5" if TalyFEM was compiled with HDF5 support
target_link_libraries(nptel_ode_ts talyfem p4est::p4est)
target_link_libraries(nptel_poisson talyfem p4est::p4est)
target_link_libraries(nptel_heat2d_ts talyfem p4est::p4est)
target_link_libraries(nptel_turing_patterns talyfem p4est::p4est)
target_link_libraries(petsc_ex45 talyfem p4est::p4est)
target_link_libraries(p4est_simple p4est::p4est)


## Try to find PETSc 3.6+
#find_package(PETSc 3.6 COMPONENTS C REQUIRED)
#
## Try to find MPI
#if (NOT DEFINED MPI_HOME)
#    set(MPI_HOME ${PETSC_DIR}/${PETSC_ARCH})
#endif()
#find_package(MPI COMPONENTS C REQUIRED)
#string(STRIP "${MPI_C_LINK_FLAGS}" MPI_C_LINK_FLAGS)
#
## List include directories
#target_include_directories(petsc_ts
#        PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>
#        PUBLIC ${PETSC_INCLUDES}
#        PUBLIC ${MPI_C_INCLUDE_PATH}
#        PUBLIC ${CONFIG++_INCLUDE_DIR}
#)
#
## List other libraries
#target_link_libraries(petsc_ts
#        LINK_PUBLIC ${PETSC_LIBRARIES} ${MPI_C_LIBRARIES} ${MPI_C_LINK_FLAGS} ${CONFIG++_LIBRARY}
#)

add_executable(heat2d src/my_examples/heat2d.cpp)
target_include_directories(heat2d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(heat2d talyfem)

add_executable(heat2d_ptc src/my_examples/heat2d_pseudo_time.cpp)
target_include_directories(heat2d_ptc PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(heat2d_ptc talyfem)

add_executable(vdp src/vanderpol/vanderpol.cpp)
target_include_directories(vdp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(vdp talyfem)

add_executable(vdp_bvp src/vanderpol/vanderpol_time_bvp.cpp)
target_include_directories(vdp_bvp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(vdp_bvp talyfem)

add_executable(vdp_bvp_ptc src/vanderpol/vanderpol_time_bvp_ptc.cpp)
target_include_directories(vdp_bvp_ptc PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(vdp_bvp_ptc talyfem)

add_executable(vdp_periodic_ptc src/vanderpol/vanderpol_time_periodic_ptc.cpp)
target_include_directories(vdp_periodic_ptc PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(vdp_periodic_ptc talyfem)

add_executable(vdp_periodic_Tunk_ptc src/vanderpol/vanderpol_time_periodic_Tunkown_ptc.cpp)
target_include_directories(vdp_periodic_Tunk_ptc PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(vdp_periodic_Tunk_ptc talyfem)

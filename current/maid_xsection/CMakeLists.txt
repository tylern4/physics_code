cmake_minimum_required(VERSION 3.5)
project(xsection C CXX Fortran)
set(VERSION 1.0.0)


set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_FLAGS "-Ofast")
set(CMAKE_C_FLAGS "-Ofast")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

add_library(FortranLib STATIC
    src/aao.F
    src/aaosig1.F
    src/cgln_amps.F
    src/daresbury.F
    src/dsigma.F
    src/helicity_amps.F
    src/interp.F
    src/legendre.F
    src/maid.F
    src/maid_lee.F
    src/maid_xsection.F
    src/multipole_amps.F
    src/read_sf_file.F
    src/splie2.F
    src/splin2.F
    src/spline.F
    src/splint.F
    src/xsection.F
    src/revinm.F
    src/fint.F
    src/lenocc.F
    src/reclog.F
    src/chklun.F
    src/dumb_functions.F
    src/chkfpk.F
    src/get_log.F
    src/get_prompt.F
    src/keep_prompt.F
    src/maidinterface.cc
    src/keep_log.F
  )



add_executable(maid_csv maid_csv.cpp)
target_link_libraries(maid_csv PUBLIC FortranLib)
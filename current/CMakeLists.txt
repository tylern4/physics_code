cmake_minimum_required(VERSION 3.2)
project(e1d)

#set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS Core Imt RIO Net Hist Graf Graf3d Gpad Tree TreePlayer Rint Postscript Matrix Physics MathCore Thread MultiProc)
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})
include(${ROOT_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/cpp)


# Materials
set(PHYS
  src/branches.cpp
  src/cuts.cpp
  src/datahandeler.cpp
  src/delta_t.cpp
  src/event.cpp
  src/fits.cpp
  src/func.cpp
  src/glob_files.cpp
  src/histogram.cpp
  src/makeHeader.cpp
  src/missing_mass.cpp
  src/photon_flux.cpp
  src/physics.cpp)
add_library(physics ${PHYS})

target_link_libraries(physics ${ROOT_LIBRARIES})

set(E1D cpp/main_cpp.cpp)
add_executable(e1d ${E1D} ${PHYS} ${ROOT_LIBRARIES})

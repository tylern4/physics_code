cmake_minimum_required(VERSION 3.8)

#find_package(PythonInterp)
#find_package(PythonLibs)
#find_package(PythonExtensions)
set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net)
include(${ROOT_USE_FILE})

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

add_library(physics STATIC
  photon_flux.cpp
  branches.cpp
  func.cpp
  datahandeler.cpp
  cuts.cpp
  skim.cpp
  histogram.cpp
  physics.cpp
  fits.cpp
  glob_files.cpp
  missing_mass.cpp
  delta_t.cpp
  reaction.cpp
  yeilds.cpp
  makeHeader.cpp
  ntuple.cpp
  fid_cuts.f
  mom_corr.cpp
  vertex_correction.f
  syncfile.cpp
  )

set(EXTRA_FLAGS "-pthread -fPIC -Wno-return-type -Wno-writable-strings -Wno-dangling-else -Wno-comment -D__ROOFIT_NOBANNER")
target_link_libraries(physics PUBLIC ${OpenMP_CXX_LIB_NAMES} ${ROOT_LIBRARIES})

target_include_directories(physics
  PRIVATE ${CMAKE_SOURCE_DIR}/include
)

install(TARGETS physics
  EXPORT PhysicsTargets
    ARCHIVE DESTINATION lib
    PUBLIC_HEADER DESTINATION .
    )


install(
    DIRECTORY
    .
    DESTINATION
    include
    FILES_MATCHING
    PATTERN "*.hpp"
)

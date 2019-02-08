cmake_minimum_required(VERSION 3.5)

find_package(PythonInterp)
find_package(PythonLibs)
find_package(PythonExtensions)
if(${PYTHON_VERSION_MAJOR} STREQUAL "2")
  set(CMAKE_CXX_STANDARD 14)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
else()
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net) #RooFit RooFitCore RooStats RootAuth
include(${ROOT_USE_FILE})

add_library(physics STATIC
  photon_flux.cpp
  branches.cpp
  func.cpp
  datahandeler.cpp
  datahandeler_mc.cpp
  cuts.cpp
  skim.cpp
  histogram.cpp
  histogram_mc.cpp
  physics.cpp
  fits.cpp
  glob_files.cpp
  missing_mass.cpp
  delta_t.cpp
  reaction.cpp
  yeilds.cpp
  )

set(EXTRA_FLAGS "-pthread -fPIC -m64 -Wno-return-type -Wno-writable-strings -Wno-dangling-else -Wno-comment -D__ROOFIT_NOBANNER")
target_link_libraries(physics ${ROOT_LIBRARIES})

target_include_directories(physics
  PRIVATE include
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

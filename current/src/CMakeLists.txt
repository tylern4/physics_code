cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net)
include(${ROOT_USE_FILE})

add_library(physics SHARED
  photon_flux.cpp
  branches.cpp
  Particle.cpp
  func.cpp
  datahandeler.cpp
  cuts.cpp
  skim.cpp
  histogram.cpp
  physics.cpp
  fits.cpp
  glob_files.cpp
  missing_mass.cpp
  Event.cpp
  delta_t.cpp
  makeHeader.cpp
  )

set(EXTRA_FLAGS "-pthread -m64 -Wno-return-type -Wno-writable-strings -Wno-dangling-else -Wno-comment")
target_link_libraries(physics ${ROOT_LIBRARIES})

target_include_directories(physics
  PRIVATE include
)

install(TARGETS physics
  EXPORT PhysicsTargets
    LIBRARY DESTINATION lib
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

cmake_minimum_required(VERSION 3.8)

# For CMake < 3.9, we need to make the target ourselves
if(NOT TARGET OpenMP::OpenMP_CXX)
    find_package(Threads REQUIRED)
    add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
    # Only works if the same flag is passed to the linker; use CMake 3.9+ otherwise (Intel, AppleClang)
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)

endif()



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
  )

set(EXTRA_FLAGS "-pthread -fPIC -Wno-return-type -Wno-writable-strings -Wno-dangling-else -Wno-comment -D__ROOFIT_NOBANNER")
target_link_libraries(physics PUBLIC OpenMP::OpenMP_CXX ${ROOT_LIBRARIES})

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

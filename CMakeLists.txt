cmake_minimum_required(VERSION 3.11)
project(DataHandler C CXX Fortran)
set(VERSION 1.0.0)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# add_subdirectory(${CMAKE_SOURCE_DIR}/external/aao_rad_lund)


list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net) #RooFit RooFitCore RooStats RootAuth
include(${ROOT_USE_FILE})

set(CMAKE_MODULE_PATH
   ${CMAKE_MODULE_PATH}
   ${PROJECT_SOURCE_DIR}/cmake
   ${PROJECT_SOURCE_DIR}/cmake/modules)

# Optimisation flags.
#
# These must be set BEFORE add_subdirectory() below, otherwise the `physics`
# library -- which is where essentially all of the analysis code lives -- is
# configured with the CMake defaults instead of these. The set() used to sit
# below add_subdirectory(), so the library silently ignored the project's
# flags while the executables (added later) picked them up.
#
# The flags are plain -O3, deliberately with no fast-math.
#
# -Ofast implies -ffinite-math-only, which lets the compiler assume no NaN or Inf
# is ever produced and so constant-folds std::isnan/std::isinf to false. This
# code uses NaN both as a "not computed yet" sentinel and as a value to be
# rejected before binning, and NaNs do reach ROOT's TH1::Fill, whose bin
# selection for NaN is unspecified. Dropping -Ofast keeps those guards.
#
# -ffp-contract=fast was tried here and reverted. It is worth nothing measurable
# (the event loop is ~94% ROOT decompression) and it changes the last digits of
# every kinematic variable: with the csv_maker e1d output, -ffp-contract=fast
# altered theta and mm2 at the ~1e-7 level in 750 of 877 rows relative to plain
# -O3. Cross sections should not depend on whether the compiler fused a
# multiply-add, so the numerics are left exactly as -O3 produces them.
#
# (With CMAKE_BUILD_TYPE=Release a later bare -O3 also overrode -Ofast on the
# command line, so the old setting was doubly inert.)
set(OPT_FLAGS "-O3")
set(CMAKE_CXX_FLAGS "${OPT_FLAGS} ${ROOT_CXX_FLAGS}")
set(CMAKE_C_FLAGS "${OPT_FLAGS} ${ROOT_CXX_FLAGS}")

# include_directories(${CMAKE_SOURCE_DIR}/external/zstr/src)
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/include)
add_subdirectory(${CMAKE_SOURCE_DIR}/src)


# include_directories(${CMAKE_SOURCE_DIR}/external/lz4xx/include)
# add_subdirectory(${CMAKE_SOURCE_DIR}/external/lz4xx)
# set_target_properties(lz4xx PROPERTIES CXX_STANDARD 11)

set(PROGRAMS
  exe/cpp/e1d.cpp
  exe/cpp/e1f.cpp
  exe/cpp/e16.cpp

  # exe/cpp/analysis.cpp
  exe/monteCarlo/e1d_mc.cpp
  exe/monteCarlo/e1f_mc.cpp
  exe/monteCarlo/analysis_mc.cpp
  exe/skim/e1d_skim.cpp
  exe/skim/e1f_skim.cpp
  exe/skim/skim_one_check_banks.cpp
  exe/csv_maker/csv_maker.cpp
  exe/csv_maker/csv_maker_PPi0.cpp
  exe/csv_maker/csv_maker_mc.cpp
  exe/ntuple_maker/ntuple_maker.cpp
  exe/golden_run/golden_run.cpp
  exe/momentumCorrections/momemtumCorrections.cpp
  exe/peakFitMomCorr/peakFitMomCorr.cpp
  exe/deltaTFitter/deltaTFitter.cpp
  exe/deltaTFitter/deltaTFitter_e1f.cpp

  # exe/fidFitter/fidFitter.cpp
  # exe/fidFitter/fidFitterPip.cpp
  511-lab/511_lab.cpp
  511-lab/make_electron_csv.cpp
)

foreach(prog ${PROGRAMS} )
    get_filename_component(progname ${prog} NAME)
    string( REPLACE ".cpp" "" progname ${progname} )
    add_executable( ${progname} ${prog} )
    set_target_properties(${progname} PROPERTIES CXX_STANDARD 17)
    target_link_libraries(${progname} PUBLIC physics ${ROOT_LIBRARIES})
    install(TARGETS ${progname} EXPORT ${PROJECT_NAME}Targets RUNTIME DESTINATION bin)
endforeach( prog ${PROGRAMS} )

#find_package(Cython)
#IF(Cython_FOUND)
#IF(0)
#  find_package(Python3 COMPONENTS NumPy REQUIRED)
  # find_package(PythonLibs 3.5 REQUIRED)
  #  set(CMAKE_CXX_STANDARD 17)
  #set(CMAKE_CXX_STANDARD_REQUIRED ON)

  #include_directories(${PYTHON_INCLUDE_DIRS} ${NumPy_INCLUDE_DIR})

  #set(EXTRA_FLAGS "${EXTRA_FLAGS} -pthread -Wno-deprecated-register")
  #set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${EXTRA_FLAGS}" )
  #include_directories(${PYTHON_INCLUDE_DIRS} ${PYTHON_NUMPY_INCLUDE_DIR} ${NumPy_INCLUDE_DIR})

  #file(GLOB CYTHON_TARGETS src/*.pyx )
  #foreach(pyx ${CYTHON_TARGETS} )
  #  get_filename_component(pyxname ${pyx} NAME)
  #  string( REPLACE ".pyx" "" pyxname ${pyxname} )
  #    add_cython_target(${pyxname} ${pyx} CXX PY3 OUTPUT_VAR _pyxname)
  #    add_library(${pyxname} MODULE ${_pyxname})
  #    python_extension_module(${pyxname})
  #    set_target_properties(${pyxname} PROPERTIES CXX_STANDARD 17)
  #    target_link_libraries(${pyxname} physics ${PYTHON_LIBRARIES} ${ROOT_LIBRARIES})
  #    install(TARGETS ${pyxname} LIBRARY DESTINATION lib)
  #endforeach( pyx ${CYTHON_TARGETS} )
  #SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  #ENDIF()


file(GLOB EXAMPLES scripts/*.cxx )
foreach(exe ${EXAMPLES} )
    get_filename_component(exename ${exe} NAME)
    string( REPLACE ".cxx" "" exename ${exename} )
    add_executable( ${exename} ${exe} )
    set_target_properties(${exename} PROPERTIES CXX_STANDARD 17)
    include_directories(${exename} PUBLIC src/include)
    target_link_libraries(${exename} PUBLIC ${ROOT_LIBRARIES})
    install(TARGETS ${exename} EXPORT ${PROJECT_NAME}Targets RUNTIME DESTINATION bin)
endforeach( exe ${EXAMPLES} )



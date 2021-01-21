cmake_minimum_required(VERSION 3.5)
project(DataHandler C CXX Fortran)
set(VERSION 1.0.0)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(${CMAKE_SOURCE_DIR}/external/aao_rad_lund)


list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net) #RooFit RooFitCore RooStats RootAuth
include(${ROOT_USE_FILE})

set(CMAKE_MODULE_PATH
   ${CMAKE_MODULE_PATH}
   ${PROJECT_SOURCE_DIR}/cmake
   ${PROJECT_SOURCE_DIR}/cmake/modules)

include_directories(${CMAKE_SOURCE_DIR}/external/zstr/src)
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/include)
add_subdirectory(${CMAKE_SOURCE_DIR}/src)


# include_directories(${CMAKE_SOURCE_DIR}/external/lz4xx/include)
# add_subdirectory(${CMAKE_SOURCE_DIR}/external/lz4xx)
# set_target_properties(lz4xx PROPERTIES CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "-Ofast ${ROOT_CXX_FLAGS}")
set(CMAKE_C_FLAGS "-Ofast ${ROOT_CXX_FLAGS}")

set(PROGRAMS
      exe/cpp/e1d.cpp
      exe/cpp/e1f.cpp
      exe/cpp/e16.cpp
      #exe/cpp/analysis.cpp
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
      exe/fidFitter/fidFitter.cpp
      )

foreach(prog ${PROGRAMS} )
    get_filename_component(progname ${prog} NAME)
    string( REPLACE ".cpp" "" progname ${progname} )
    add_executable( ${progname} ${prog} )
    set_target_properties(${progname} PROPERTIES CXX_STANDARD 17)
    target_link_libraries(${progname} PUBLIC physics ${ROOT_LIBRARIES})
    install(TARGETS ${progname} EXPORT ${PROJECT_NAME}Targets RUNTIME DESTINATION bin)
endforeach( prog ${PROGRAMS} )

find_package(Cython)
IF(Cython_FOUND)
  find_package(PythonInterp 3.5 REQUIRED)
  find_package(PythonLibs 3.5 REQUIRED)
  find_package(PythonExtensions REQUIRED)
  find_package(NumPy REQUIRED)
  set(CMAKE_CXX_STANDARD 14)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)

  include_directories(${PYTHON_INCLUDE_DIRS} ${NumPy_INCLUDE_DIR})

  set(EXTRA_FLAGS "${EXTRA_FLAGS} -pthread -Wno-deprecated-register")
  set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${EXTRA_FLAGS}" )
  include_directories(${PYTHON_INCLUDE_DIRS} ${PYTHON_NUMPY_INCLUDE_DIR} ${NumPy_INCLUDE_DIR})

  file(GLOB CYTHON_TARGETS src/*.pyx )
  foreach(pyx ${CYTHON_TARGETS} )
    get_filename_component(pyxname ${pyx} NAME)
    string( REPLACE ".pyx" "" pyxname ${pyxname} )
      add_cython_target(${pyxname} ${pyx} CXX PY3 OUTPUT_VAR _pyxname)
      add_library(${pyxname} MODULE ${_pyxname})
      python_extension_module(${pyxname})
      set_target_properties(${pyxname} PROPERTIES CXX_STANDARD 17)
      target_link_libraries(${pyxname} physics ${PYTHON_LIBRARIES} ${ROOT_LIBRARIES})
      install(TARGETS ${pyxname} LIBRARY DESTINATION lib)
  endforeach( pyx ${CYTHON_TARGETS} )
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
ENDIF()


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



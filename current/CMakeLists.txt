cmake_minimum_required(VERSION 3.5)
project(DataHandeler C CXX Fortran)
set(VERSION 1.0.0)

if(UNIX AND NOT APPLE)
    define_property(-DLINUX)
endif()

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net) #RooFit RooFitCore RooStats RootAuth
include(${ROOT_USE_FILE})

set(CMAKE_MODULE_PATH
   ${CMAKE_MODULE_PATH}
   ${PROJECT_SOURCE_DIR}/cmake
   ${PROJECT_SOURCE_DIR}/cmake/modules)

include_directories("src")
add_subdirectory("src")

set(CMAKE_CXX_FLAGS "-Ofast ${ROOT_CXX_FLAGS}")
set(CMAKE_C_FLAGS "-Ofast ${ROOT_CXX_FLAGS}")

set(PROGRAMS 
      cpp/e1d.cpp 
      cpp/e1f.cpp 
      cpp/analysis.cpp
      511-lab/511_lab.cpp 
      monteCarlo/e1d_mc.cpp 
      monteCarlo/e1f_mc.cpp 
      monteCarlo/analysis_mc.cpp 
      skim/e1d_skim.cpp 
      skim/e1f_skim.cpp 
      csv_maker/csv_maker.cpp 
      csv_maker/csv_maker_mc.cpp 
      ntuple_maker/ntuple_maker.cpp
      golden_run/golden_run.cpp
      momentumCorrections/momemtumCorrections.cpp
      peakFitMomCorr/peakFitMomCorr.cpp
      deltaTFitter/deltaTFitter.cpp
      )
    
foreach(prog ${PROGRAMS} )
    get_filename_component(progname ${prog} NAME)
    string( REPLACE ".cpp" "" progname ${progname} )
    add_executable( ${progname} ${prog} )
    target_link_libraries(${progname} PUBLIC physics ${ROOT_LIBRARIES})
    install(TARGETS ${progname} EXPORT ${PROJECT_NAME}Targets RUNTIME DESTINATION bin)
endforeach( prog ${PROGRAMS} )

find_package(Cython)
IF(Cython_FOUND)
  find_package(PythonInterp REQUIRED)
  find_package(PythonLibs REQUIRED)
  find_package(PythonExtensions REQUIRED)
  find_package(NumPy REQUIRED)
  if(${PYTHON_VERSION_MAJOR} STREQUAL "2")
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
  else()
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
  endif()

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
    include_directories(${exename} PUBLIC src/include)
    target_link_libraries(${exename} PUBLIC ${ROOT_LIBRARIES})
    install(TARGETS ${exename} EXPORT ${PROJECT_NAME}Targets RUNTIME DESTINATION bin)
endforeach( exe ${EXAMPLES} )
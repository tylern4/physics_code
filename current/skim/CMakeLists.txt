cmake_minimum_required(VERSION 3.6)
project(hipoio VERSION 0.0.1 DESCRIPTION "skim")

set(CMAKE_MODULE_PATH
   ${CMAKE_MODULE_PATH}
   ${PROJECT_SOURCE_DIR}/cmake
   ${PROJECT_SOURCE_DIR}/cmake/modules)

set(HEADERS
  ../src/skim.hpp
  )


include(cmake/os.cmake)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT REQUIRED COMPONENTS Tree)
include(${ROOT_USE_FILE})
include_directories("../src" ${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIR})
include(cmake/root.cmake)
set(lib_HEADERS ${HEADERS})

find_package(PythonInterp)
find_package(PythonLibs)
find_package(PythonExtensions)
find_package(Cython REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_cython_target(skim skim.pyx CXX PY3 OUTPUT_VAR _skim)
include_directories(${PYTHON_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}" )

add_library(skim MODULE ${_skim})
python_extension_module(skim)

#target_link_libraries(${PYTHON_LIBRARIES})

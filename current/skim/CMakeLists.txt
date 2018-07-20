cmake_minimum_required(VERSION 3.5)
project(pySkim)
set(CMAKE_MODULE_PATH
   ${CMAKE_MODULE_PATH}
   ${PROJECT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}/etc/cmake)
include(cmake/os.cmake)
include(cmake/root.cmake)

find_package(PythonInterp)
find_package(PythonLibs)
find_package(PythonExtensions)
find_package(Cython REQUIRED)
find_package(ROOT REQUIRED COMPONENTS Tree)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(phySkim STATIC skim.cpp)


add_cython_target(skim skim.pyx CXX PY2 OUTPUT_VAR _skim)
include_directories(_skim PUBLIC
  ${PROJECT_SOURCE_DIR}
  ${PYTHON_INCLUDE_DIRS}
  ${ROOT_INCLUDE_DIRS}
  ../src)

add_library(skim MODULE ${_skim})
add_dependencies(skim phySkim)
target_link_libraries(skim phySkim ${PYTHON_LIBRARIES} ${ROOT_LIBRARIES})
python_extension_module(skim)

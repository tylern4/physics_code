cmake_minimum_required(VERSION 3.5)
project(DataHandeler)
set(VERSION 0.0.1)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net)
include(${ROOT_USE_FILE})

set(CMAKE_MODULE_PATH
   ${CMAKE_MODULE_PATH}
   ${PROJECT_SOURCE_DIR}/cmake
   ${PROJECT_SOURCE_DIR}/cmake/modules)

include_directories("current/src")
add_subdirectory("current/src")

find_package(PythonInterp)
find_package(PythonLibs)
find_package(PythonExtensions)
find_package(Cython REQUIRED)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_cython_target(handeler handeler.pyx CXX PY2 OUTPUT_VAR _handeler)
include_directories(${PYTHON_INCLUDE_DIRS})


set(EXTRA_FLAGS "-pthread -m64 -Wno-deprecated-register")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${EXTRA_FLAGS}" )

add_library(handeler MODULE ${_handeler})
python_extension_module(handeler)

target_link_libraries(handeler physics ${PYTHON_LIBRARIES} ${ROOT_LIBRARIES})

install(TARGETS handeler
        LIBRARY DESTINATION lib)

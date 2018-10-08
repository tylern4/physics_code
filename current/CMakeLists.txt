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

include_directories("src")
add_subdirectory("src")

# e1d
add_executable(e1d cpp/main_cpp.cpp)
target_link_libraries(e1d physics ${ROOT_LIBRARIES})
install(TARGETS e1d
  RUNTIME DESTINATION bin)

#Skim
add_executable(e1d_skim skim/main_skim.cpp)
target_link_libraries(e1d_skim physics ${ROOT_LIBRARIES})
install(TARGETS e1d_skim
  RUNTIME DESTINATION bin)


find_package(PythonInterp)
find_package(PythonLibs)
find_package(PythonExtensions)
find_package(Cython)

IF(Cython_FOUND)
  if(${PYTHON_VERSION_MAJOR} STREQUAL "2")
    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
  else()
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
  endif()
  set(EXTRA_FLAGS "-pthread -m64 -Wno-deprecated-register")
  set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${EXTRA_FLAGS}" )
  include_directories(${PYTHON_INCLUDE_DIRS} ${PYTHON_NUMPY_INCLUDE_DIR} ${NumPy_INCLUDE_DIR})

  # Python DataHandeler
  add_cython_target(handeler src/handeler.pyx CXX PY3 OUTPUT_VAR _handeler)
  add_library(handeler MODULE ${_handeler})
  python_extension_module(handeler)
  target_link_libraries(handeler physics ${PYTHON_LIBRARIES} ${ROOT_LIBRARIES})
  install(TARGETS handeler
          LIBRARY DESTINATION lib)

  # Python Skim
  add_cython_target(skim src/skim.pyx CXX PY3 OUTPUT_VAR _skim)
  add_library(skim MODULE ${_skim})
  python_extension_module(skim)
  target_link_libraries(skim physics ${PYTHON_LIBRARIES} ${ROOT_LIBRARIES})
  install(TARGETS skim
          LIBRARY DESTINATION lib)

  # Python h10
  add_cython_target(h10 src/h10.pyx CXX PY3 OUTPUT_VAR _h10)
  add_library(h10 MODULE ${_h10})
  python_extension_module(h10)
  target_link_libraries(h10 physics ${PYTHON_LIBRARIES} ${ROOT_LIBRARIES})
  install(TARGETS h10
          LIBRARY DESTINATION lib)
ENDIF()

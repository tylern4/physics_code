cmake_minimum_required(VERSION 3.5)
project(DataHandeler)
set(VERSION 0.0.1)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net) #RooFit RooFitCore RooStats RootAuth
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

# lab_511
add_executable(lab_511 511-lab/511_lab.cpp)
target_link_libraries(lab_511 physics ${ROOT_LIBRARIES})
install(TARGETS lab_511
  RUNTIME DESTINATION bin)

# e1d_mc
add_executable(e1d_mc monteCarlo/main_mc.cpp)
target_link_libraries(e1d_mc physics ${ROOT_LIBRARIES})
install(TARGETS e1d_mc
    RUNTIME DESTINATION bin)

#Skim
add_executable(e1d_skim skim/main_skim.cpp)
target_link_libraries(e1d_skim physics ${ROOT_LIBRARIES})
install(TARGETS e1d_skim
  RUNTIME DESTINATION bin)

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

  set(EXTRA_FLAGS "-pthread -Wno-deprecated-register -D__ROOFIT_NOBANNER")
  set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${EXTRA_FLAGS}" )
  include_directories(${PYTHON_INCLUDE_DIRS} ${PYTHON_NUMPY_INCLUDE_DIR} ${NumPy_INCLUDE_DIR})

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

  # Python h10
  add_cython_target(yeilds src/yeilds.pyx CXX PY3 OUTPUT_VAR _yeilds)
  add_library(yeilds MODULE ${_yeilds})
  python_extension_module(yeilds)
  target_link_libraries(yeilds physics ${PYTHON_LIBRARIES} ${ROOT_LIBRARIES})
  install(TARGETS yeilds LIBRARY DESTINATION lib)
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
ENDIF()

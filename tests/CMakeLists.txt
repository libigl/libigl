cmake_minimum_required(VERSION 3.1)
project(libigl_tutorials)
message(STATUS "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER}")
message(STATUS "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")

### conditionally compile certain modules depending on libraries found on the system
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake)
#find_package(CGAL QUIET COMPONENTS Core)
find_package(MATLAB QUIET)
find_package(MOSEK)

### libIGL options: choose between header only and compiled static library
option(LIBIGL_USE_STATIC_LIBRARY "Use LibIGL as static library" ON)
option(LIBIGL_WITH_EMBREE      "Use Embree"         OFF)

### libIGL options: choose your dependencies (by default everything is OFF, in this example we need the viewer)
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" ON)
option(LIBIGL_WITH_CGAL              "Use CGAL"           OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"         OFF)
option(LIBIGL_WITH_CORK              "Use Cork"           OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"         OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK             "Use MOSEK"          OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         OFF)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           OFF)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          OFF)
option(LIBIGL_WITH_PNG               "Use PNG"            OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"       OFF)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  OFF)
option(LIBIGL_WITH_XML               "Use XML"            OFF)
option(LIBIGL_WITH_PYTHON            "Use Python"         OFF)
### End

### Adding libIGL: choose the path to your local copy libIGL
include(libigl)

### Download data
igl_download_test_data()
set(IGL_TEST_DATA ${LIBIGL_EXTERNAL}/../tests/data)

### Download Google unit test framework.
igl_download_googletest()

SET(TEST_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR}) 
INCLUDE_DIRECTORIES(${TEST_ROOT_DIR})

# Set TEST_DIR definition
ADD_DEFINITIONS(-DLIBIGL_DATA_DIR="${IGL_TEST_DATA}")

# Add googletest googlemock support
ENABLE_TESTING()
ADD_SUBDIRECTORY(
  ${LIBIGL_EXTERNAL}/googletest/googlemock
  ${CMAKE_CURRENT_BINARY_DIR}/gtest)
SET(GTEST_BOTH_LIBRARIES gtest gmock)
INCLUDE_DIRECTORIES(${gmock_SOURCE_DIR})
INCLUDE_DIRECTORIES(${gmock_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${gtest_SOURCE_DIR})
INCLUDE_DIRECTORIES(${gtest_SOURCE_DIR}/include)

# Process code in each subdirectories: add in decreasing order of complexity
# (last added will run first and those should be the fastest tests)
IF(LIBIGL_WITH_MOSEK)
  ADD_SUBDIRECTORY(${TEST_ROOT_DIR}/include/igl/mosek)
ENDIF()

IF(LIBIGL_WITH_CGAL)
  ADD_SUBDIRECTORY(${TEST_ROOT_DIR}/include/igl/copyleft/boolean)
  ADD_SUBDIRECTORY(${TEST_ROOT_DIR}/include/igl/copyleft/cgal)
ENDIF()

IF(LIBIGL_WITH_TETGEN)
  ADD_SUBDIRECTORY(${TEST_ROOT_DIR}/include/igl/copyleft/tetgen)
ENDIF()

IF(LIBIGL_WITH_COMISO)
  ADD_SUBDIRECTORY(${TEST_ROOT_DIR}/include/igl/copyleft/comiso)
ENDIF()

ADD_SUBDIRECTORY(${TEST_ROOT_DIR}/include/igl)

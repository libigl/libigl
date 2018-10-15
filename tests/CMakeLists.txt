cmake_minimum_required(VERSION 3.1)
project(libigl_tests)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake)		
		
### Adding libIGL: choose the path to your local copy libIGL		
if(NOT TARGET igl_common)		
  include(libigl)		
else()		
  include(LibiglDownloadExternal)		
endif()

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

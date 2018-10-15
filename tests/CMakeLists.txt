SET(TEST_ROOT_DIR ${PROJECT_SOURCE_DIR}/tests) 
INCLUDE_DIRECTORIES(${TEST_ROOT_DIR})

# Set TEST_DIR definition
ADD_DEFINITIONS(-DLIBIGL_DATA_DIR="${LIBIGL_DATA}")

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

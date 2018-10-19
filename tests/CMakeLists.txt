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
igl_download_catch2()

SET(TEST_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})
list(APPEND CMAKE_MODULE_PATH ${LIBIGL_EXTERNAL}/Catch2/contrib)
INCLUDE_DIRECTORIES(${TEST_ROOT_DIR})




# Add catch2
add_library(catch INTERFACE)
target_include_directories(catch SYSTEM INTERFACE ${LIBIGL_EXTERNAL}/catch2/single_include)

# Create test executable
add_executable(libigl_tests main.cpp test_common.h)
target_link_libraries(libigl_tests PUBLIC igl::core catch)
# Set TEST_DIR definition
set(DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data/")
target_compile_definitions(libigl_tests PUBLIC -DLIBIGL_DATA_DIR="${IGL_TEST_DATA}")


# Process code in each subdirectories: add in decreasing order of complexity
# (last added will run first and those should be the fastest tests)
IF(LIBIGL_WITH_MOSEK)
  file(GLOB TEST_SRC_FILES ${TEST_ROOT_DIR}/include/igl/mosek/*.cpp)
  file(GLOB TEST_INC_FILES ${TEST_ROOT_DIR}/include/igl/mosek/*.h ${TEST_ROOT_DIR}/include/igl/mosek/*.inl)
  target_sources(libigl_tests PRIVATE ${TEST_SRC_FILES} ${TEST_INC_FILES})

  target_link_libraries(libigl_tests PUBLIC igl::mosek)
ENDIF()

IF(LIBIGL_WITH_CGAL)
  file(GLOB TEST_SRC_FILES ${TEST_ROOT_DIR}/include/igl/copyleft/boolean/*.cpp ${TEST_ROOT_DIR}/include/igl/copyleft/cgal/*.cpp)
  file(GLOB TEST_INC_FILES ${TEST_ROOT_DIR}/include/igl/copyleft/boolean/*.h ${TEST_ROOT_DIR}/include/igl/copyleft/cgal/*.h ${TEST_ROOT_DIR}/include/igl/copyleft/boolean/*.inl ${TEST_ROOT_DIR}/include/igl/copyleft/cgal/*.inl)
  target_sources(libigl_tests PRIVATE ${TEST_SRC_FILES} ${TEST_INC_FILES})

  target_link_libraries(libigl_tests PUBLIC igl::igl::cgal)
ENDIF()

IF(LIBIGL_WITH_TETGEN)
  file(GLOB TEST_SRC_FILES ${TEST_ROOT_DIR}/include/igl/copyleft/tetgen/*.cpp)
  file(GLOB TEST_INC_FILES ${TEST_ROOT_DIR}/include/igl/copyleft/tetgen/*.h ${TEST_ROOT_DIR}/include/igl/copyleft/tetgen/*.inl)
  target_sources(libigl_tests PRIVATE ${TEST_SRC_FILES} ${TEST_INC_FILES})

  target_link_libraries(libigl_tests PUBLIC igl::tetgen)
ENDIF()

IF(LIBIGL_WITH_COMISO)
  file(GLOB TEST_SRC_FILES ${TEST_ROOT_DIR}/include/igl/copyleft/comiso/*.cpp)
  file(GLOB TEST_INC_FILES ${TEST_ROOT_DIR}/include/igl/copyleft/comiso/*.h ${TEST_ROOT_DIR}/include/igl/copyleft/comiso/*.inl)
  target_sources(libigl_tests PRIVATE ${TEST_SRC_FILES} ${TEST_INC_FILES})

  target_link_libraries(libigl_tests PUBLIC igl::comiso)
ENDIF()


file(GLOB TEST_SRC_FILES ${TEST_ROOT_DIR}/include/igl/*.cpp)
file(GLOB TEST_INC_FILES ${TEST_ROOT_DIR}/include/igl/*.h ${TEST_ROOT_DIR}/include/igl/*.inl)
target_sources(libigl_tests PRIVATE ${TEST_SRC_FILES} ${TEST_INC_FILES})




# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
include(Catch)
catch_discover_tests(libigl_tests)







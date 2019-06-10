cmake_minimum_required(VERSION 3.1)
project(libigl_tests)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake)

### Adding libIGL: choose the path to your local copy libIGL
if(NOT TARGET igl::core)
  include(libigl)
else()
  include(LibiglDownloadExternal)
endif()

### Download data
igl_download_test_data()
set(IGL_TEST_DATA ${LIBIGL_EXTERNAL}/../tests/data)

### Download Catch2 unit test framework
igl_download_catch2()
list(APPEND CMAKE_MODULE_PATH ${LIBIGL_EXTERNAL}/catch2/contrib)


# Add catch2
add_subdirectory(${LIBIGL_EXTERNAL}/catch2 catch2)


# Create test executable
add_executable(libigl_tests main.cpp test_common.h)
target_link_libraries(libigl_tests PUBLIC igl::core Catch2::Catch2)
target_include_directories(libigl_tests PUBLIC ${CMAKE_CURRENT_LIST_DIR})

# Set DATA_DIR definition
set(DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data/")
target_compile_definitions(libigl_tests PUBLIC -DLIBIGL_DATA_DIR="${IGL_TEST_DATA}")


# Process code in each subdirectories: add in decreasing order of complexity
# (last added will run first and those should be the fastest tests)
if(LIBIGL_WITH_MOSEK)
  file(GLOB TEST_SRC_FILES ./include/igl/mosek/*.cpp)
  file(GLOB TEST_INC_FILES ./include/igl/mosek/*.h ./include/igl/mosek/*.inl)
  target_sources(libigl_tests PRIVATE ${TEST_SRC_FILES} ${TEST_INC_FILES})

  target_link_libraries(libigl_tests PUBLIC igl::mosek)
endif()

if(LIBIGL_WITH_CGAL)
  file(GLOB TEST_SRC_FILES ./include/igl/copyleft/boolean/*.cpp ./include/igl/copyleft/cgal/*.cpp)
  file(GLOB TEST_INC_FILES ./include/igl/copyleft/boolean/*.h ./include/igl/copyleft/cgal/*.h ./include/igl/copyleft/boolean/*.inl ./include/igl/copyleft/cgal/*.inl)
  target_sources(libigl_tests PRIVATE ${TEST_SRC_FILES} ${TEST_INC_FILES})

  target_link_libraries(libigl_tests PUBLIC igl::cgal)
  igl_copy_cgal_dll(libigl_tests)
endif()

if(LIBIGL_WITH_TETGEN)
  file(GLOB TEST_SRC_FILES ./include/igl/copyleft/tetgen/*.cpp)
  file(GLOB TEST_INC_FILES ./include/igl/copyleft/tetgen/*.h ./include/igl/copyleft/tetgen/*.inl)
  target_sources(libigl_tests PRIVATE ${TEST_SRC_FILES} ${TEST_INC_FILES})

  target_link_libraries(libigl_tests PUBLIC igl::tetgen)
endif()

if(LIBIGL_WITH_COMISO)
  file(GLOB TEST_SRC_FILES ./include/igl/copyleft/comiso/*.cpp)
  file(GLOB TEST_INC_FILES ./include/igl/copyleft/comiso/*.h ./include/igl/copyleft/comiso/*.inl)
  target_sources(libigl_tests PRIVATE ${TEST_SRC_FILES} ${TEST_INC_FILES})

  target_link_libraries(libigl_tests PUBLIC igl::comiso)
endif()

if(LIBIGL_WITH_EMBREE)
  file(GLOB TEST_SRC_FILES ./include/igl/embree/*.cpp)
  file(GLOB TEST_INC_FILES ./include/igl/embree/*.h ./include/igl/embree/*.inl)
  target_sources(libigl_tests PRIVATE ${TEST_SRC_FILES} ${TEST_INC_FILES})

  target_link_libraries(libigl_tests PUBLIC igl::embree)
endif()

if(LIBIGL_WITH_PREDICATES)
  file(GLOB TEST_SRC_FILES ./include/igl/predicates/*.cpp)
  file(GLOB TEST_INC_FILES ./include/igl/predicates/*.h ./include/igl/predicates/*.inl)
  target_sources(libigl_tests PRIVATE ${TEST_SRC_FILES} ${TEST_INC_FILES})

  target_link_libraries(libigl_tests PUBLIC igl::predicates igl::triangle)
endif()

file(GLOB TEST_SRC_FILES ./include/igl/*.cpp)
file(GLOB TEST_INC_FILES ./include/igl/*.h ./include/igl/*.inl)
target_sources(libigl_tests PRIVATE ${TEST_SRC_FILES} ${TEST_INC_FILES})


# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
include(Catch)
catch_discover_tests(libigl_tests)

cmake_minimum_required(VERSION 3.16)
project(use-libigl)

include(CTest)
include(FetchContent)
include("./UseLibigl.cmake")

find_package(Eigen3 CONFIG REQUIRED)
find_package(libigl CONFIG REQUIRED
  COMPONENTS core copyleft_cgal
  NO_DEFAULT_PATH
  PATHS "${libigl_ROOT}")

check_subproject_dirs(libigl_tests_data)
check_subproject_dirs(Catch2)
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY "${Catch2_SOURCE_DIR}"
  GIT_TAG "HEAD")

FetchContent_MakeAvailable(Catch2)

# Create tests that use the libigl library found above.

create_use_libigl_test(
  IGL_TARGET igl::core
  IGL_TEST_SOURCES cat.cpp)

create_use_libigl_test(
  IGL_TARGET igl_copyleft::cgal
  IGL_TEST_SOURCES copyleft/cgal/triangulate.cpp)

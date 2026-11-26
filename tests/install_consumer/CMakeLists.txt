# Minimal downstream project that consumes an *installed* libigl through find_package().
#
# CI builds and installs libigl with one compiler, then builds this project with another (e.g.
# MSVC-built libigl consumed by clang++). It fails if the installed targets carry usage
# requirements that only libigl's own compiler accepts (see libigl#2590).
#
# Configure with -Dlibigl_DIR=<prefix>/lib/cmake/igl.
cmake_minimum_required(VERSION 3.16)
project(libigl_install_consumer CXX)

# libigl installs Eigen's headers and EigenTargets.cmake but no Eigen3Config.cmake; this directory
# provides a FindEigen3.cmake that loads them, for libigl-config.cmake's find_dependency(Eigen3).
list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

find_package(libigl CONFIG REQUIRED)

add_executable(consumer main.cpp)
target_link_libraries(consumer PRIVATE igl::core)

enable_testing()
add_test(NAME consumer COMMAND consumer)

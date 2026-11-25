cmake_minimum_required(VERSION 3.20)
project(exact_point_predicates LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
find_package(Threads REQUIRED)
get_filename_component(LIBIGL_POINT_SOURCE "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
foreach(required PREDICATES_SOURCE_DIR EIGEN_INCLUDE_DIR CATCH2_INCLUDE_DIR)
  if(NOT DEFINED ${required})
    message(FATAL_ERROR "Set ${required} to an existing checkout/include directory")
  endif()
endforeach()
set(FETCHCONTENT_SOURCE_DIR_PREDICATES "${PREDICATES_SOURCE_DIR}" CACHE PATH "" FORCE)
set(FETCHCONTENT_SOURCE_DIR_EIGEN "${EIGEN_INCLUDE_DIR}" CACHE PATH "" FORCE)
set(LIBIGL_PREDICATES ON CACHE BOOL "" FORCE)
add_subdirectory("${LIBIGL_POINT_SOURCE}" libigl)
add_executable(test_exact_points
  "${LIBIGL_POINT_SOURCE}/tests/main.cpp"
  "${LIBIGL_POINT_SOURCE}/tests/include/igl/predicates/exact_point_triangle.cpp"
  "${LIBIGL_POINT_SOURCE}/tests/include/igl/predicates/exact_point_mesh.cpp"
  "${LIBIGL_POINT_SOURCE}/tests/include/igl/predicates/exact_ray_triangle.cpp")
target_include_directories(test_exact_points PRIVATE "${LIBIGL_POINT_SOURCE}/tests" "${CATCH2_INCLUDE_DIR}")
target_compile_definitions(test_exact_points PRIVATE LIBIGL_DATA_DIR="${CMAKE_CURRENT_BINARY_DIR}/data")
target_link_libraries(test_exact_points PRIVATE igl::predicates)
enable_testing()
add_test(NAME exact_points COMMAND test_exact_points)
set_tests_properties(exact_points PROPERTIES TIMEOUT 120 ENVIRONMENT "IGL_NUM_THREADS=4")
add_library(exact_point_static STATIC
  "${LIBIGL_POINT_SOURCE}/include/igl/default_num_threads.cpp"
  "${LIBIGL_POINT_SOURCE}/include/igl/predicates/exactinit.cpp"
  "${LIBIGL_POINT_SOURCE}/include/igl/predicates/exact_point_triangle_closest_feature.cpp"
  "${LIBIGL_POINT_SOURCE}/include/igl/predicates/exact_point_triangle_squared_distance_compare.cpp"
  "${LIBIGL_POINT_SOURCE}/include/igl/predicates/exact_point_mesh_closest_triangle.cpp"
  "${LIBIGL_POINT_SOURCE}/include/igl/predicates/exact_point_mesh_classify.cpp"
  "${LIBIGL_POINT_SOURCE}/include/igl/predicates/exact_mesh_feature_exposure.cpp"
  "${LIBIGL_POINT_SOURCE}/include/igl/predicates/point_mesh_closest_point.cpp"
  "${LIBIGL_POINT_SOURCE}/include/igl/predicates/approximate_signed_distance_exact_sign.cpp")
target_compile_definitions(exact_point_static PUBLIC IGL_STATIC_LIBRARY)
target_include_directories(exact_point_static PUBLIC "${LIBIGL_POINT_SOURCE}/include" "${EIGEN_INCLUDE_DIR}")
target_link_libraries(exact_point_static PUBLIC predicates Threads::Threads)
if(MSVC)
  target_compile_options(exact_point_static PUBLIC /fp:strict)
else()
  target_compile_options(exact_point_static PUBLIC -fno-fast-math -ffp-contract=off)
endif()
add_executable(static_smoke static_smoke.cpp)
target_link_libraries(static_smoke PRIVATE exact_point_static)
add_test(NAME exact_point_static_api COMMAND static_smoke)
set_tests_properties(exact_point_static_api PROPERTIES TIMEOUT 30)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_executable(audit_allocations allocations.cpp)
  target_link_libraries(audit_allocations PRIVATE igl::predicates)
  target_link_options(audit_allocations PRIVATE -Wl,--wrap=predicates_expansion_sum -Wl,--wrap=nextafter)
  add_test(NAME exact_point_allocations COMMAND audit_allocations)
  set_tests_properties(exact_point_allocations PROPERTIES TIMEOUT 30)
endif()
add_executable(benchmark_exact_points benchmark.cpp)
target_link_libraries(benchmark_exact_points PRIVATE igl::predicates)
add_executable(profile_exact_point_search profile_search.cpp)
target_link_libraries(profile_exact_point_search PRIVATE igl::predicates)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_executable(audit_exact_point_search profile_search.cpp)
  target_compile_definitions(audit_exact_point_search PRIVATE PROFILE_EXPANSIONS)
  target_link_libraries(audit_exact_point_search PRIVATE igl::predicates)
  target_link_options(audit_exact_point_search PRIVATE -Wl,--wrap=predicates_expansion_sum)
endif()

# Build the headless tutorial against existing dependencies and data.
if(DEFINED TUTORIAL_DATA_DIR)
  add_executable(1008_SignedDistanceProfile
    "${LIBIGL_POINT_SOURCE}/tutorial/1008_SignedDistanceProfile/main.cpp")
  target_compile_definitions(1008_SignedDistanceProfile PRIVATE TUTORIAL_SHARED_PATH="${TUTORIAL_DATA_DIR}")
  target_link_libraries(1008_SignedDistanceProfile PRIVATE igl::predicates)
endif()

add_executable(profile_point_filters profile_filters.cpp)
target_link_libraries(profile_point_filters PRIVATE igl::predicates)
add_test(NAME predicates_patch
  COMMAND ${CMAKE_COMMAND}
    "-DPREDICATES_SOURCE_DIR=${PREDICATES_SOURCE_DIR}"
    "-DPATCH_DIR=${LIBIGL_POINT_SOURCE}/cmake/recipes/external"
    "-DTEST_DIR=${CMAKE_CURRENT_BINARY_DIR}/predicates-patch-test"
    -P "${LIBIGL_POINT_SOURCE}/tests/cmake/predicates_patch.cmake")
set_tests_properties(predicates_patch PROPERTIES TIMEOUT 60)

cmake_minimum_required(VERSION 3.20)
project(exact_ray_comparison LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
get_filename_component(LIBIGL_EXPERIMENT_SOURCE "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
foreach(required BASELINE_SOURCE_DIR PREDICATES_SOURCE_DIR EIGEN_INCLUDE_DIR CATCH2_INCLUDE_DIR)
  if(NOT DEFINED ${required})
    message(FATAL_ERROR "Set ${required} to an existing checkout/include directory")
  endif()
endforeach()
set(FETCHCONTENT_SOURCE_DIR_PREDICATES "${PREDICATES_SOURCE_DIR}" CACHE PATH "" FORCE)
set(FETCHCONTENT_SOURCE_DIR_EIGEN "${EIGEN_INCLUDE_DIR}" CACHE PATH "" FORCE)
set(LIBIGL_PREDICATES ON CACHE BOOL "" FORCE)
add_subdirectory("${LIBIGL_EXPERIMENT_SOURCE}" libigl)
add_library(baseline_adapter STATIC adapter.cpp)
target_compile_definitions(baseline_adapter PRIVATE EXACT_RAYS_BASELINE=1)
target_include_directories(baseline_adapter PRIVATE "${BASELINE_SOURCE_DIR}/include" "${EIGEN_INCLUDE_DIR}")
add_library(shewchuk_adapter STATIC adapter.cpp)
target_include_directories(shewchuk_adapter PRIVATE "${LIBIGL_EXPERIMENT_SOURCE}/include" "${EIGEN_INCLUDE_DIR}")
target_link_libraries(shewchuk_adapter PRIVATE igl::predicates)
add_executable(compare_exact_rays compare.cpp)
target_include_directories(compare_exact_rays PRIVATE "${EIGEN_INCLUDE_DIR}")
target_link_libraries(compare_exact_rays PRIVATE baseline_adapter shewchuk_adapter)
add_executable(audit_allocations allocations.cpp)
target_include_directories(audit_allocations PRIVATE "${EIGEN_INCLUDE_DIR}")
target_link_libraries(audit_allocations PRIVATE baseline_adapter shewchuk_adapter)

add_executable(test_exact_rays
  "${LIBIGL_EXPERIMENT_SOURCE}/tests/main.cpp"
  "${LIBIGL_EXPERIMENT_SOURCE}/tests/include/igl/predicates/predicates.cpp"
  "${LIBIGL_EXPERIMENT_SOURCE}/tests/include/igl/predicates/orient3d.cpp"
  "${LIBIGL_EXPERIMENT_SOURCE}/tests/include/igl/predicates/exact_ray_mesh_intersect.cpp"
  "${LIBIGL_EXPERIMENT_SOURCE}/tests/include/igl/predicates/exact_ray_triangle.cpp")
target_include_directories(test_exact_rays PRIVATE
  "${LIBIGL_EXPERIMENT_SOURCE}/include" "${LIBIGL_EXPERIMENT_SOURCE}/tests"
  "${EIGEN_INCLUDE_DIR}" "${CATCH2_INCLUDE_DIR}")
target_compile_definitions(test_exact_rays PRIVATE LIBIGL_DATA_DIR="${CMAKE_CURRENT_BINARY_DIR}/data")
target_link_libraries(test_exact_rays PRIVATE igl::predicates)

# Exercise the actual explicit instantiations without building unrelated modules.
add_library(exact_ray_static STATIC
  "${LIBIGL_EXPERIMENT_SOURCE}/include/igl/predicates/exactinit.cpp"
  "${LIBIGL_EXPERIMENT_SOURCE}/include/igl/predicates/exact_ray_triangle_intersect.cpp"
  "${LIBIGL_EXPERIMENT_SOURCE}/include/igl/predicates/exact_ray_triangle_hit_compare.cpp"
  "${LIBIGL_EXPERIMENT_SOURCE}/include/igl/predicates/exact_ray_mesh_intersect.cpp"
  "${LIBIGL_EXPERIMENT_SOURCE}/include/igl/eytzinger_aabb.cpp"
  "${LIBIGL_EXPERIMENT_SOURCE}/include/igl/eytzinger_aabb_ray_intersection.cpp"
  "${LIBIGL_EXPERIMENT_SOURCE}/include/igl/median.cpp"
  "${LIBIGL_EXPERIMENT_SOURCE}/include/igl/matrix_to_list.cpp")
target_compile_definitions(exact_ray_static PUBLIC IGL_STATIC_LIBRARY)
target_include_directories(exact_ray_static PUBLIC "${LIBIGL_EXPERIMENT_SOURCE}/include" "${EIGEN_INCLUDE_DIR}")
target_link_libraries(exact_ray_static PUBLIC predicates)
add_executable(static_smoke static_smoke.cpp)
target_link_libraries(static_smoke PRIVATE exact_ray_static)
add_executable(core_smoke core_smoke.cpp)
target_link_libraries(core_smoke PRIVATE igl::core)

foreach(target baseline_adapter shewchuk_adapter compare_exact_rays test_exact_rays exact_ray_static)
  if(MSVC)
    target_compile_options(${target} PRIVATE /fp:strict)
  else()
    target_compile_options(${target} PRIVATE -fno-fast-math -ffp-contract=off)
  endif()
endforeach()
enable_testing()
add_test(NAME predicate_regression COMMAND test_exact_rays "~exact_ray_mesh_intersect: agrees_with_double")
add_test(NAME baseline_comparison COMMAND compare_exact_rays --accuracy)
add_test(NAME expansion_allocations COMMAND audit_allocations)
add_test(NAME static_api COMMAND static_smoke)
add_test(NAME core_without_predicates COMMAND core_smoke)

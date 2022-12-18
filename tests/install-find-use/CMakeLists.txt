set(libigl_ROOT "${CMAKE_CURRENT_BINARY_DIR}/install")

FetchContent_GetProperties(
  Catch2
  SOURCE_DIR Catch2_SOURCE_DIR
  BINARY_DIR Catch2_BINARY_DIR)

add_test(
  NAME install-libigl
  COMMAND
  "${CMAKE_COMMAND}"
  --install "${PROJECT_BINARY_DIR}"
  --prefix "${libigl_ROOT}"
  --config $<CONFIG>)

set_tests_properties(install-libigl PROPERTIES FIXTURES_SETUP install-fixture)

function(create_find_use_test igl_target)
  string(REPLACE "::" "_" igl_target "${igl_target}")
  set(test_name "find-use-${igl_target}")

  add_test(
    NAME ${test_name}
    COMMAND
    "${CMAKE_CTEST_COMMAND}"
    --verbose
    --output-on-failure
    --build-noclean
    --build-generator "${CMAKE_GENERATOR}"
    --build-config $<CONFIG>
    --build-target "use-${igl_target}"
    --build-and-test
      "${CMAKE_CURRENT_SOURCE_DIR}/use-libigl"
      "${CMAKE_CURRENT_BINARY_DIR}/use-libigl"
    --build-options
      "-Dlibigl_ROOT:PATH=${libigl_ROOT}"
      "-DLIBIGL_DATA_DIR:PATH=${LIBIGL_DATA_DIR}"
      "-DEigen3_ROOT:PATH=${Eigen3_ROOT}"
      "-DCatch2_SOURCE_DIR:PATH=${Catch2_SOURCE_DIR}"
      "-Dlibigl_tests_data_SOURCE_DIR:PATH=${libigl_tests_data_SOURCE_DIR}"
      "-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}"
    --test-command "${CMAKE_CTEST_COMMAND}" --verbose --output-on-failure)

  set_tests_properties(${test_name} PROPERTIES
    RESOURCE_LOCK "use-libigl" # opposed to separate build directories
    FIXTURES_REQUIRED install-fixture)
endfunction()


create_find_use_test(igl::core)

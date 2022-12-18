include(CMakeParseArguments)

function(check_subproject_dirs project)
  if(NOT ${project}_SOURCE_DIR)
    message(FATAL_ERROR
      "'${project}_SOURCE_DIR' must be specified. It should be populated by FetchContent when the "
      "specific version of ${project} is acquired by the libigl build. Received "
      "${project}_SOURCE_DIR='${${project}_SOURCE_DIR}'")
  endif()
endfunction()

function(create_use_libigl_test)
  cmake_parse_arguments(ARGS
    ""                            # options
    IGL_TARGET                    # one-value
    IGL_TEST_SOURCES;TEST_SOURCES # multi-value
    "${ARGN}")

  if(NOT DEFINED ARGS_IGL_TARGET)
    message(FATAL_ERROR "The argument 'IGL_TARGET' must be provided to ${CMAKE_CURRENT_FUNCTION}")
  endif()


  set(igl_tests_root "${PROJECT_SOURCE_DIR}/../..")
  set(test_exec "use-${ARGS_IGL_TARGET}")
  string(REPLACE "::" "_" test_exec "${test_exec}")

  # to steal existing igl test code
  set(igl_sources "${ARGS_IGL_TEST_SOURCES}")
  list(TRANSFORM igl_sources PREPEND "${igl_tests_root}/include/igl/")

  set(test_sources "${ARGS_TEST_SOURCES}")
  list(TRANSFORM test_sources PREPEND "${PROJECT_SOURCE_DIR}")

  add_executable(${test_exec}
    "${igl_sources}"
    "${test_sources}"
    "${igl_tests_root}/main.cpp")

  target_compile_definitions(${test_exec} PRIVATE LIBIGL_DATA_DIR="${libigl_tests_data_SOURCE_DIR}")
  target_include_directories(${test_exec} PRIVATE "${igl_tests_root}")
  target_link_libraries(${test_exec} PRIVATE
    Eigen3::Eigen
    Catch2::Catch2
    ${ARGS_IGL_TARGET})

  add_test(NAME ${test_exec} COMMAND ${test_exec})
endfunction()

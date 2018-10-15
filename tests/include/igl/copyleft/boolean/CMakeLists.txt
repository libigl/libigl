set(TEST_SRC_FILES main.cpp mesh_boolean.cpp)

add_executable(igl_boolean_tests ${TEST_SRC_FILES} ${TEST_INC_FILES})
target_link_libraries(igl_boolean_tests igl::core igl::cgal gtest_main)
add_test(NAME run_igl_boolean_tests COMMAND igl_boolean_tests)

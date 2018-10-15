set(TEST_SRC_FILES bbw.cpp main.cpp)

add_executable(igl_mosek_tests ${TEST_SRC_FILES} ${TEST_INC_FILES})
target_link_libraries(igl_mosek_tests igl::core igl::mosek gtest_main)
add_test(NAME run_igl_mosek_tests COMMAND igl_mosek_tests)

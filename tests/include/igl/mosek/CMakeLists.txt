file(GLOB TEST_SRC_FILES *.cpp main.cpp)
file(GLOB TEST_INC_FILES *.h *.inl)

add_executable(igl_mosek_tests ${TEST_SRC_FILES} ${TEST_INC_FILES})
target_link_libraries(igl_mosek_tests igl::core igl::mosek gtest_main)
add_test(NAME run_igl_mosek_tests COMMAND igl_mosek_tests)

file(GLOB TEST_SRC_FILES *.cpp main.cpp)
file(GLOB TEST_INC_FILES *.h *.inl)

add_executable(igl_tetgen_tests ${TEST_SRC_FILES} ${TEST_INC_FILES})
target_link_libraries(igl_tetgen_tests igl::core igl::tetgen gtest_main)
add_test(NAME run_igl_tetgen_tests COMMAND igl_tetgen_tests)

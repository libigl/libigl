# Check if CGAL is available
file(GLOB TEST_SRC_FILES *.cpp main.cpp)
file(GLOB TEST_INC_FILES *.h *.inl)

add_executable(igl_boolean_tests ${TEST_SRC_FILES} ${TEST_INC_FILES})
target_link_libraries(igl_boolean_tests igl::core igl::cgal gtest_main)
add_test(NAME run_igl_boolean_tests COMMAND igl_boolean_tests)

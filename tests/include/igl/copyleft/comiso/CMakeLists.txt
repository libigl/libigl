file(GLOB TEST_SRC_FILES *.cpp main.cpp)
file(GLOB TEST_INC_FILES *.h *.inl)

add_executable(igl_comiso_tests ${TEST_SRC_FILES} ${TEST_INC_FILES})
target_link_libraries(igl_comiso_tests igl::core igl::comiso gtest_main)
add_test(NAME run_igl_comiso_tests COMMAND igl_comiso_tests)

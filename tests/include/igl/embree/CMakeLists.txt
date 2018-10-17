file(GLOB TEST_SRC_FILES *.cpp main.cpp)
file(GLOB TEST_INC_FILES *.h *.inl)

add_executable(igl_embree_tests ${TEST_SRC_FILES} ${TEST_INC_FILES})
target_link_libraries(igl_embree_tests igl::core igl::embree gtest_main)
add_test(NAME run_igl_embree_tests COMMAND igl_embree_tests)

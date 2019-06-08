file(GLOB TEST_SRC_FILES *.cpp main.cpp)
file(GLOB TEST_INC_FILES *.h *.inl)

add_executable(igl_tests ${TEST_SRC_FILES} ${TEST_INC_FILES})
target_link_libraries(igl_tests igl::core gtest_main)
add_test(NAME run_igl_tests COMMAND igl_tests)

if(MSVC)
  # Silencing some compile warnings
  target_compile_options(igl_tests PRIVATE 
    # Type conversion warnings. These can be fixed with some effort and possibly more verbose code.
    /wd4267 # conversion from 'size_t' to 'type', possible loss of data
    /wd4244 # conversion from 'type1' to 'type2', possible loss of data
    /wd4018 # signed/unsigned mismatch
    /wd4305 # truncation from 'double' to 'float'
  )
endif(MSVC)
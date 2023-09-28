if(TARGET igl::tests_data)
    return()
endif()

message(STATUS "Third-party: creating target 'igl::tests_data'")

include(FetchContent)
FetchContent_Declare(
    libigl_tests_data
    GIT_REPOSITORY https://github.com/libigl/libigl-tests-data
    GIT_TAG        19cedf96d70702d8b3a83eb27934780c542356fe
)
FetchContent_MakeAvailable(libigl_tests_data)

add_library(igl_tests_data INTERFACE)
add_library(igl::tests_data ALIAS igl_tests_data)

target_compile_definitions(igl_tests_data INTERFACE LIBIGL_DATA_DIR=\"${libigl_tests_data_SOURCE_DIR}\")

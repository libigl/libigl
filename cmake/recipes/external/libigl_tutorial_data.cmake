if(TARGET igl::tutorial_data)
    return()
endif()

message(STATUS "Third-party: creating target 'igl::tutorial_data'")

include(FetchContent)
FetchContent_Declare(
    libigl_tutorial_data
    GIT_REPOSITORY https://github.com/libigl/libigl-tutorial-data
    GIT_TAG        644dd4104843b6d736745d9dafbd70bf8d175648
)
FetchContent_MakeAvailable(libigl_tutorial_data)

add_library(igl_tutorial_data INTERFACE)
add_library(igl::tutorial_data ALIAS igl_tutorial_data)

target_compile_definitions(igl_tutorial_data INTERFACE "-DTUTORIAL_SHARED_PATH=\"${libigl_tutorial_data_SOURCE_DIR}\"")

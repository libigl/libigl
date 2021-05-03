if(TARGET igl::tutorial_data)
    return()
endif()

message(STATUS "Third-party: creating target 'igl::tutorial_data'")

include(FetchContent)
FetchContent_Declare(
    libigl_tutorial_tata
    GIT_REPOSITORY https://github.com/libigl/libigl-tutorial-data
    GIT_TAG        c1f9ede366d02e3531ecbaec5e3769312f31cccd
)
FetchContent_MakeAvailable(libigl_tutorial_tata)

add_library(igl_tutorial_data INTERFACE)
add_library(igl::tutorial_data ALIAS igl_tutorial_data)

target_compile_definitions(igl_tutorial_data INTERFACE "-DTUTORIAL_SHARED_PATH=\"${libigl_tutorial_tata_SOURCE_DIR}\"")

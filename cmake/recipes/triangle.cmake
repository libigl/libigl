if(TARGET triangle::triangle)
    return()
endif()

message(STATUS "Third-party: creating target 'triangle::triangle'")

include(FetchContent)
FetchContent_Declare(
    triangle
    GIT_REPOSITORY https://github.com/libigl/triangle.git
    GIT_TAG        5a70326574b34d6a51d9eaf6a9f78813657ee108
)

FetchContent_MakeAvailable(triangle)
add_library(triangle::triangle ALIAS triangle)

target_include_directories(triangle INTERFACE "${triangle_SOURCE_DIR}")

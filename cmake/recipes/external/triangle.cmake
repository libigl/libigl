if(TARGET triangle::triangle)
    return()
endif()

message(STATUS "Third-party: creating target 'triangle::triangle'")

include(FetchContent)
FetchContent_Declare(
    triangle
    GIT_REPOSITORY https://github.com/libigl/triangle.git
    GIT_TAG        62f02db9ab4ff4b62d5ff82a77c8ea458c84c23a
)

FetchContent_MakeAvailable(triangle)
add_library(triangle::triangle ALIAS triangle)

target_include_directories(triangle INTERFACE "${triangle_SOURCE_DIR}")

set_target_properties(triangle PROPERTIES FOLDER ThirdParty)

if(TARGET triangle::triangle)
    return()
endif()

message(STATUS "Third-party: creating target 'triangle::triangle'")

include(FetchContent)
FetchContent_Declare(
    triangle
    GIT_REPOSITORY https://github.com/libigl/triangle.git
    GIT_TAG        efe1f1b672a3885266c42feaba94a8bd73db450c
)

FetchContent_MakeAvailable(triangle)
add_library(triangle::triangle ALIAS triangle)

target_include_directories(triangle INTERFACE "${triangle_SOURCE_DIR}")

set_target_properties(triangle PROPERTIES FOLDER ThirdParty)

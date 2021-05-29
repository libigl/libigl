if(TARGET triangle::triangle)
    return()
endif()

message(STATUS "Third-party: creating target 'triangle::triangle'")

include(FetchContent)
FetchContent_Declare(
    triangle
    GIT_REPOSITORY https://github.com/libigl/triangle.git
    GIT_TAG        3ee6cac2230f0fe1413879574f741c7b6da11221
)

FetchContent_MakeAvailable(triangle)
add_library(triangle::triangle ALIAS triangle)

target_include_directories(triangle INTERFACE "${triangle_SOURCE_DIR}")

set_target_properties(triangle PROPERTIES FOLDER ThirdParty)

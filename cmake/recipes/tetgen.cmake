if(TARGET tetgen::tetgen)
    return()
endif()

message(STATUS "Third-party: creating target 'tetgen::tetgen'")

include(FetchContent)
FetchContent_Declare(
    tetgen
    GIT_REPOSITORY https://github.com/libigl/tetgen.git
    GIT_TAG        f6804b31625cf88e7e1b76fc5c6d6441716d3fc1
)

FetchContent_MakeAvailable(tetgen)
add_library(tetgen::tetgen ALIAS tetgen)

target_include_directories(tetgen INTERFACE "${tetgen_SOURCE_DIR}")

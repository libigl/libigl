if(TARGET tetgen::tetgen)
    return()
endif()

message(STATUS "Third-party: creating target 'tetgen::tetgen'")

include(FetchContent)
FetchContent_Declare(
    tetgen
    GIT_REPOSITORY https://github.com/libigl/tetgen.git
    GIT_TAG        4f3bfba3997f20aa1f96cfaff604313a8c2c85b6
)

FetchContent_MakeAvailable(tetgen)
add_library(tetgen::tetgen ALIAS tetgen)

target_include_directories(tetgen INTERFACE "${tetgen_SOURCE_DIR}")

set_target_properties(tetgen PROPERTIES FOLDER ThirdParty)

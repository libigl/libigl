if(TARGET tetgen::tetgen)
    return()
endif()

message(STATUS "Third-party: creating target 'tetgen::tetgen'")

include(FetchContent)
FetchContent_Declare(
    tetgen
    GIT_REPOSITORY https://github.com/libigl/tetgen.git
    GIT_TAG        e05aca7df74e3f531bc35733ed87d36d437266c5
)

FetchContent_MakeAvailable(tetgen)
add_library(tetgen::tetgen ALIAS tetgen)

target_include_directories(tetgen INTERFACE "${tetgen_SOURCE_DIR}")

set_target_properties(tetgen PROPERTIES FOLDER ThirdParty)

if(TARGET glad::glad)
    return()
endif()

message(STATUS "Third-party: creating target 'glad::glad'")

include(FetchContent)
FetchContent_Declare(
    glad
    GIT_REPOSITORY https://github.com/KeithBallard/libigl-glad.git
    GIT_TAG        09a93ab
)

FetchContent_MakeAvailable(glad)
add_library(glad::glad ALIAS glad)

set_target_properties(glad PROPERTIES FOLDER ThirdParty)

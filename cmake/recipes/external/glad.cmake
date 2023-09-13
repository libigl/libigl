if(TARGET glad::glad)
    return()
endif()

message(STATUS "Third-party: creating target 'glad::glad'")

include(FetchContent)
FetchContent_Declare(
    glad
    GIT_REPOSITORY https://github.com/libigl/libigl-glad.git
    GIT_TAG        b9490d66dcc94ef0aeaded0bebd0ed8c515d4404
)

FetchContent_MakeAvailable(glad)
add_library(glad::glad ALIAS glad)

set_target_properties(glad PROPERTIES FOLDER ThirdParty)

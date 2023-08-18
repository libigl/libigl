if(TARGET glad::glad)
    return()
endif()

message(STATUS "Third-party: creating target 'glad::glad'")

include(FetchContent)
FetchContent_Declare(
    glad
    GIT_REPOSITORY https://github.com/libigl/libigl-glad.git
    GIT_TAG        ceef55fcd08bdd16e985370a99cfb60e69623221
)

FetchContent_MakeAvailable(glad)
add_library(glad::glad ALIAS glad)

set_target_properties(glad PROPERTIES FOLDER ThirdParty)

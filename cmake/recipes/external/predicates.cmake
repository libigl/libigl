if(TARGET predicates::predicates)
    return()
endif()

message(STATUS "Third-party: creating target 'predicates::predicates'")

include(FetchContent)
FetchContent_Declare(
    predicates
    GIT_REPOSITORY https://github.com/libigl/libigl-predicates.git
    GIT_TAG        488242fa2b1f98a9c5bd1441297fb4a99a6a9ae4
)

FetchContent_MakeAvailable(predicates)
add_library(predicates::predicates ALIAS predicates)

set_target_properties(predicates PROPERTIES FOLDER ThirdParty)

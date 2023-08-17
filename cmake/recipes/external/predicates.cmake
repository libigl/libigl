if(TARGET predicates::predicates)
    return()
endif()

message(STATUS "Third-party: creating target 'predicates::predicates'")

include(FetchContent)
FetchContent_Declare(
    predicates
    GIT_REPOSITORY https://github.com/KeithBallard/libigl-predicates.git
    GIT_TAG        d6c04da
)

FetchContent_MakeAvailable(predicates)
add_library(predicates::predicates ALIAS predicates)

set_target_properties(predicates PROPERTIES FOLDER ThirdParty)

if(TARGET predicates::predicates)
    return()
endif()

message(STATUS "Third-party: creating target 'predicates::predicates'")

include(FetchContent)
FetchContent_Declare(
    predicates
    GIT_REPOSITORY https://github.com/libigl/libigl-predicates.git
    GIT_TAG        decb7bc1260e689cbe008109e3cc5d3a5a433aea
)

FetchContent_MakeAvailable(predicates)
add_library(predicates::predicates ALIAS predicates)

set_target_properties(predicates PROPERTIES FOLDER ThirdParty)

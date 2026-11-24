if(TARGET predicates::predicates)
    return()
endif()

message(STATUS "Third-party: creating target 'predicates::predicates'")

include(FetchContent)
FetchContent_Declare(
    predicates
    GIT_REPOSITORY https://github.com/libigl/libigl-predicates.git
    GIT_TAG        decb7bc1260e689cbe008109e3cc5d3a5a433aea
    # Experimental public wrappers for Shewchuk's existing expansion arithmetic.
    PATCH_COMMAND ${CMAKE_COMMAND} -DSOURCE_DIR=<SOURCE_DIR>
        -P ${CMAKE_CURRENT_LIST_DIR}/predicates-patch.cmake
)

FetchContent_MakeAvailable(predicates)
add_library(predicates::predicates ALIAS predicates)

set_target_properties(predicates PROPERTIES FOLDER ThirdParty)

if(TARGET cork::cork)
    return()
endif()

message(STATUS "Third-party: creating target 'cork::cork'")

include(FetchContent)
FetchContent_Declare(
    cork
    GIT_REPOSITORY https://github.com/libigl/cork.git
    GIT_TAG        27ad8a285838f5a480d856429e39d3d56d4338f9
)

include(gmp)

FetchContent_MakeAvailable(cork)
add_library(cork::cork ALIAS cork)

target_include_directories(cork INTERFACE "${cork_SOURCE_DIR}/src")

set_target_properties(cork PROPERTIES FOLDER ThirdParty)
set_target_properties(cork-bin PROPERTIES FOLDER ThirdParty)

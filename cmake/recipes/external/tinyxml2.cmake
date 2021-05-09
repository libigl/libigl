if(TARGET tinyxml2::tinyxml2)
    return()
endif()

message(STATUS "Third-party: creating target 'tinyxml2::tinyxml2'")

# TODO: Update tinyxml2 version + use FetchContent_MakeAvailable with upstream CMake code
include(FetchContent)
FetchContent_Declare(
    tinyxml2
    GIT_REPOSITORY https://github.com/leethomason/tinyxml2.git
    GIT_TAG        d175e9de0be0d4db75d0a8cf065599a435a87eb6
)
FetchContent_GetProperties(tinyxml2)
if(NOT tinyxml2_POPULATED)
    FetchContent_Populate(tinyxml2)
endif()

add_library(tinyxml2 STATIC ${tinyxml2_SOURCE_DIR}/tinyxml2.cpp ${tinyxml2_SOURCE_DIR}/tinyxml2.h)
add_library(tinyxml2::tinyxml2 ALIAS tinyxml2)
target_include_directories(tinyxml2 PUBLIC ${tinyxml2_SOURCE_DIR})
set_target_properties(tinyxml2 PROPERTIES DEFINE_SYMBOL "TINYXML2_EXPORT")

set_target_properties(tinyxml2 PROPERTIES FOLDER ThirdParty)

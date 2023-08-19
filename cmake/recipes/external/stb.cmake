if(TARGET stb::stb)
    return()
endif()

message(STATUS "Third-party: creating target 'stb::stb'")

include(FetchContent)
FetchContent_Declare(
    stb
    GIT_REPOSITORY https://github.com/nothings/stb.git
    GIT_TAG f67165c2bb2af3060ecae7d20d6f731173485ad0
)
FetchContent_MakeAvailable(stb)

# Generate implementation file
file(WRITE "${stb_BINARY_DIR}/stb_image.cpp.in" [[
    #define STB_IMAGE_IMPLEMENTATION
    #include <stb_image.h>

    #define STB_IMAGE_WRITE_IMPLEMENTATION
    #include <stb_image_write.h>
]])

configure_file(${stb_BINARY_DIR}/stb_image.cpp.in ${stb_BINARY_DIR}/stb_image.cpp COPYONLY)

# Define stb library
add_library(stb ${stb_BINARY_DIR}/stb_image.cpp)
add_library(stb::stb ALIAS stb)

target_include_directories(stb PUBLIC "${stb_SOURCE_DIR}")

set_target_properties(stb PROPERTIES FOLDER ThirdParty)

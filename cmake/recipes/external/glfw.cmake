if(TARGET glfw::glfw)
    return()
endif()

message(STATUS "Third-party: creating target 'glfw::glfw'")

include(FetchContent)
FetchContent_Declare(
    glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG 3327050ca66ad34426a82c217c2d60ced61526b7
)

option(GLFW_BUILD_EXAMPLES "Build the GLFW example programs" OFF)
option(GLFW_BUILD_TESTS "Build the GLFW test programs" OFF)
option(GLFW_BUILD_DOCS "Build the GLFW documentation" OFF)
option(GLFW_INSTALL "Generate installation target" OFF)
option(GLFW_VULKAN_STATIC "Use the Vulkan loader statically linked into application" OFF)
FetchContent_MakeAvailable(glfw)

add_library(glfw::glfw ALIAS glfw)

set_target_properties(glfw PROPERTIES FOLDER ThirdParty)

# Warning config
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    target_compile_options(glfw PRIVATE
        "-Wno-missing-field-initializers"
        "-Wno-objc-multiple-method-names"
    )
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    target_compile_options(glfw PRIVATE
        "-Wno-missing-field-initializers"
        "-Wno-objc-multiple-method-names"
    )
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options(glfw PRIVATE
        "-Wno-missing-field-initializers"
        "-Wno-sign-compare"
        "-Wno-unused-parameter"
    )
endif()

if(TARGET imgui_fonts::imgui_fonts)
    return()
endif()

message(STATUS "Third-party: creating target 'imgui_fonts::imgui_fonts'")

include(FetchContent)
FetchContent_Declare(
    imgui_fonts
    GIT_REPOSITORY https://github.com/libigl/libigl-imgui.git
    GIT_TAG        7e1053e750b0f4c129b046f4e455243cb7f804f3
)
FetchContent_GetProperties(imgui_fonts)
if(NOT imgui_fonts_POPULATED)
    FetchContent_Populate(imgui_fonts)
endif()

add_library(imgui_fonts INTERFACE)
add_library(imgui_fonts::imgui_fonts ALIAS imgui_fonts)

include(GNUInstallDirs)
target_include_directories(imgui_fonts SYSTEM INTERFACE
    $<BUILD_INTERFACE:${imgui_fonts_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

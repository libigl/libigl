if(TARGET igl::imgui_fonts)
    return()
endif()

message(STATUS "Third-party: creating target 'igl::imgui_fonts'")

include(FetchContent)
FetchContent_Declare(
    libigl_imgui_fonts
    GIT_REPOSITORY https://github.com/libigl/libigl-imgui.git
    GIT_TAG        7e1053e750b0f4c129b046f4e455243cb7f804f3
)
FetchContent_GetProperties(libigl_imgui_fonts)
if(NOT libigl_imgui_fonts_POPULATED)
    FetchContent_Populate(libigl_imgui_fonts)
endif()

add_library(igl_imgui_fonts INTERFACE)
add_library(igl::imgui_fonts ALIAS igl_imgui_fonts)

include(GNUInstallDirs)
target_include_directories(igl_imgui_fonts SYSTEM INTERFACE
    $<BUILD_INTERFACE:${libigl_imgui_fonts_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

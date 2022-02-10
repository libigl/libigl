# 1. Define module
igl_add_library(igl_imgui)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_imgui ${IGL_SCOPE}
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# 3. Target sources
file(GLOB INC_FILES "${PROJECT_SOURCE_DIR}/include/igl/opengl/glfw/imgui/*.h")
file(GLOB SRC_FILES "${PROJECT_SOURCE_DIR}/include/igl/opengl/glfw/imgui/*.cpp")
igl_target_sources(igl_imgui ${INC_FILES} ${SRC_FILES})

# 4. Dependencies
include(imgui)
include(imguizmo)
include(libigl_imgui_fonts)
igl_include(glfw)
target_link_libraries(igl_imgui ${IGL_SCOPE}
    igl::core
    igl::glfw
    imgui::imgui
    imguizmo::imguizmo
    igl::imgui_fonts
)

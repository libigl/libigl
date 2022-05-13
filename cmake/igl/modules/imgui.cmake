# 1. Define module
igl_add_library(igl_imgui)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_imgui ${IGL_SCOPE}
  $<BUILD_INTERFACE:${libigl_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )

# 3. Target sources

igl_glob_sources("${libigl_SOURCE_DIR}/include/igl/opengl/glfw/imgui/" SRC_FILES)
igl_target_sources(igl_imgui ${SRC_FILES})

# 4. Install target & headers
igl_install(igl_imgui ${SRC_FILES})

# 5. Dependencies
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

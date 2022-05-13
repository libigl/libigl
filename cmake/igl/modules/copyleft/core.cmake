# 1. Define module
igl_add_library(igl_copyleft_core)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_copyleft_core ${IGL_SCOPE}
  $<BUILD_INTERFACE:${libigl_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )

# 3. Target sources
igl_glob_sources("${libigl_SOURCE_DIR}/include/igl/copyleft/" SRC_FILES)
igl_target_sources(igl_copyleft_core ${SRC_FILES})

# 4. Install target & headers
igl_install(igl_copyleft_core ${SRC_FILES})

# 5. Dependencies
target_link_libraries(igl_copyleft_core ${IGL_SCOPE}
  igl::core
  )

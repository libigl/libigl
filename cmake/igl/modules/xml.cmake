# 1. Define module
igl_add_library(igl_xml)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_xml ${IGL_SCOPE}
  $<BUILD_INTERFACE:${libigl_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )

# 3. Target sources
igl_glob_sources("${libigl_SOURCE_DIR}/include/igl/xml/" SRC_FILES)
igl_target_sources(igl_xml ${SRC_FILES})

# 4. Install target & headers
igl_install(igl_xml ${SRC_FILES})

# 5. Dependencies
include(tinyxml2)
target_link_libraries(igl_xml ${IGL_SCOPE}
  igl::core
  tinyxml2::tinyxml2
  )

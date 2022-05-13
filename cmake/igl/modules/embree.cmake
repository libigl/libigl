# 1. Define module
igl_add_library(igl_embree)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_embree ${IGL_SCOPE}
  $<BUILD_INTERFACE:${libigl_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )

# 3. Target sources
igl_glob_sources("${libigl_SOURCE_DIR}/include/igl/embree/" SRC_FILES)
igl_target_sources(igl_embree ${SRC_FILES})

# 4. Install target & headers
igl_install(igl_embree ${SRC_FILES})

# 5. Dependencies
include(embree)
target_link_libraries(igl_embree ${IGL_SCOPE}
  igl::core
  embree::embree
  )

# 6. Unit tests
file(GLOB SRC_FILES "${libigl_SOURCE_DIR}/tests/include/igl/embree/*.cpp")
igl_add_test(igl_embree ${SRC_FILES})

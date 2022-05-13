# 1. Define module
igl_add_library(igl_copyleft_tetgen)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_copyleft_tetgen ${IGL_SCOPE}
  $<BUILD_INTERFACE:${libigl_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )

# 3. Target sources
igl_glob_sources("${libigl_SOURCE_DIR}/include/igl/copyleft/tetgen/" SRC_FILES)
igl_target_sources(igl_copyleft_tetgen ${SRC_FILES})

# 4. Install target & headers
igl_install(igl_copyleft_tetgen ${SRC_FILES})

# 5. Dependencies
include(tetgen)
igl_include(copyleft core)
target_link_libraries(igl_copyleft_tetgen ${IGL_SCOPE}
  igl::core
  igl_copyleft::core
  tetgen::tetgen
  )

# 6. Unit tests
file(GLOB SRC_FILES "${libigl_SOURCE_DIR}/tests/include/igl/copyleft/tetgen/*.cpp")
igl_add_test(igl_copyleft_tetgen ${SRC_FILES})

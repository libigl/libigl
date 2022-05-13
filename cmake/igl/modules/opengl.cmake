# 1. Define module
igl_add_library(igl_opengl)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_opengl ${IGL_SCOPE}
  $<BUILD_INTERFACE:${libigl_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )

# 3. Target sources
igl_glob_sources("${libigl_SOURCE_DIR}/include/igl/opengl/" SRC_FILES)
igl_target_sources(igl_opengl ${SRC_FILES})

# 4. Install target & headers
igl_install(igl_opengl ${SRC_FILES})

# 5. Dependencies
include(glad)
find_package(OpenGL REQUIRED OPTIONAL_COMPONENTS OpenGL)
target_link_libraries(igl_opengl ${IGL_SCOPE}
  igl::core
  glad::glad
  # Link against OpenGL::OpenGL if available, or fallback to OpenGL::GL
  $<IF:$<TARGET_EXISTS:OpenGL::OpenGL>,OpenGL::OpenGL,OpenGL::GL>
  )

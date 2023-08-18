# 1. Define module
igl_add_library(igl_stb)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_stb ${IGL_SCOPE}
    $<BUILD_INTERFACE:${libigl_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# 3. Target sources
file(GLOB INC_FILES "${libigl_SOURCE_DIR}/include/igl/stb/*.h")
file(GLOB SRC_FILES "${libigl_SOURCE_DIR}/include/igl/stb/*.cpp")
if(LIBIGL_OPENGL)
  message(STATUS "Including igl/opengl/stb support")
  file(GLOB OPENGL_INC_FILES "${libigl_SOURCE_DIR}/include/igl/opengl/stb/*.h")
  file(GLOB OPENGL_SRC_FILES "${libigl_SOURCE_DIR}/include/igl/opengl/stb/*.cpp")
  list(APPEND INC_FILES ${OPENGL_INC_FILES})
  list(APPEND SRC_FILES ${OPENGL_SRC_FILES})
endif()
igl_target_sources(igl_stb ${INC_FILES} ${SRC_FILES})

# 4. Dependencies
include(stb)
target_link_libraries(igl_stb ${IGL_SCOPE}
    igl::core
    stb::stb
)

if(LIBIGL_OPENGL)
  igl_include(opengl)
  target_link_libraries(igl_stb ${IGL_SCOPE}
      igl::opengl
  )
endif()

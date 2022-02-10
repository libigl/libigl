# 1. Define module
igl_add_library(igl_png)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_png ${IGL_SCOPE}
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# 3. Target sources
file(GLOB INC_FILES "${PROJECT_SOURCE_DIR}/include/igl/png/*.h")
file(GLOB SRC_FILES "${PROJECT_SOURCE_DIR}/include/igl/png/*.cpp")
igl_target_sources(igl_png ${INC_FILES} ${SRC_FILES})

# 4. Dependencies
include(stb)
igl_include(opengl)
target_link_libraries(igl_png ${IGL_SCOPE}
    igl::core
    igl::opengl
    stb::stb
)

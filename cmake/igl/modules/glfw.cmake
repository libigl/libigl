# 1. Define module
igl_add_library(igl_glfw)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_glfw ${IGL_SCOPE}
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# 3. Target sources
file(GLOB INC_FILES "${PROJECT_SOURCE_DIR}/include/igl/opengl/glfw/*.h")
file(GLOB SRC_FILES "${PROJECT_SOURCE_DIR}/include/igl/opengl/glfw/*.cpp")
igl_target_sources(igl_glfw ${INC_FILES} ${SRC_FILES})

# 4. Dependencies
include(glfw)
igl_include(opengl)
target_link_libraries(igl_glfw ${IGL_SCOPE}
    igl::core
    igl::opengl
    glfw::glfw
)

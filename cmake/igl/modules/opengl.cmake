# 1. Define module
igl_add_library(igl_opengl)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_opengl ${IGL_SCOPE}
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# 3. Target sources
file(GLOB INC_FILES "${PROJECT_SOURCE_DIR}/include/igl/opengl/*.h")
file(GLOB SRC_FILES "${PROJECT_SOURCE_DIR}/include/igl/opengl/*.cpp")
igl_target_sources(igl_opengl ${INC_FILES} ${SRC_FILES})

# 4. Dependencies
include(glad)
find_package(OpenGL REQUIRED OPTIONAL_COMPONENTS OpenGL)
target_link_libraries(igl_opengl ${IGL_SCOPE}
    igl::core
    glad::glad
    # Link against OpenGL::OpenGL if available, or fallback to OpenGL::GL
    $<IF:$<TARGET_EXISTS:OpenGL::OpenGL>,OpenGL::OpenGL,OpenGL::GL>
)

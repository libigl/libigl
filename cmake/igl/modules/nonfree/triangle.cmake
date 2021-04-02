# 1. Define module
igl_add_library(igl_nonfree_triangle)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_nonfree_triangle PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# 3. Target sources
file(GLOB INC_FILES "${PROJECT_SOURCE_DIR}/include/igl/triangle/*.h")
file(GLOB SRC_FILES "${PROJECT_SOURCE_DIR}/include/igl/triangle/*.cpp")
igl_target_sources(igl_nonfree_triangle ${INC_FILES} ${SRC_FILES})

# 4. Dependencies
include(triangle)
target_link_libraries(igl_nonfree_triangle ${IGL_SCOPE}
    igl::core
    triangle::triangle
)

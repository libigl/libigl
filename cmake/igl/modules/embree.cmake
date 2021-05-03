# 1. Define module
igl_add_library(igl_embree)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_embree ${IGL_SCOPE}
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# 3. Target sources
file(GLOB INC_FILES "${PROJECT_SOURCE_DIR}/include/igl/embree/*.h")
file(GLOB SRC_FILES "${PROJECT_SOURCE_DIR}/include/igl/embree/*.cpp")
igl_target_sources(igl_embree ${INC_FILES} ${SRC_FILES})

# 4. Dependencies
include(embree)
target_link_libraries(igl_embree ${IGL_SCOPE}
    igl::core
    embree::embree
)

# 5. Unit tests
file(GLOB SRC_FILES "${PROJECT_SOURCE_DIR}/tests/include/igl/embree/*.cpp")
igl_add_test(igl_embree ${SRC_FILES})

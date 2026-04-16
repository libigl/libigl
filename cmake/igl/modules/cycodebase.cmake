# 1. Define module
igl_add_library(igl_cycodebase)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_cycodebase ${IGL_SCOPE}
    $<BUILD_INTERFACE:${libigl_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# 3. Target sources
file(GLOB INC_FILES "${libigl_SOURCE_DIR}/include/igl/cycodebase/*.h")
file(GLOB SRC_FILES "${libigl_SOURCE_DIR}/include/igl/cycodebase/*.cpp")
igl_target_sources(igl_cycodebase ${INC_FILES} ${SRC_FILES})

# 4. Dependencies
include(cycodebase)
target_link_libraries(igl_cycodebase ${IGL_SCOPE}
    igl::core
    cyCodeBase::cyCodeBase
)

# 5. Unit tests
file(GLOB SRC_FILES "${libigl_SOURCE_DIR}/tests/include/igl/cycodebase/*.cpp")
igl_add_test(igl_cycodebase ${SRC_FILES})


# 1. Define module
igl_add_library(igl_copyleft_core)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_copyleft_core PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# 3. Target sources
file(GLOB INC_FILES "${PROJECT_SOURCE_DIR}/include/igl/copyleft/*.h")
file(GLOB SRC_FILES "${PROJECT_SOURCE_DIR}/include/igl/copyleft/*.cpp")
igl_target_sources(igl_copyleft_core ${INC_FILES} ${SRC_FILES})

# 4. Dependencies
target_link_libraries(igl_copyleft_core ${IGL_SCOPE}
    igl::core
)

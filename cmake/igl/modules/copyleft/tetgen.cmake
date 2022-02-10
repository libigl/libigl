# 1. Define module
igl_add_library(igl_copyleft_tetgen)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_copyleft_tetgen ${IGL_SCOPE}
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# 3. Target sources
file(GLOB INC_FILES "${PROJECT_SOURCE_DIR}/include/igl/copyleft/tetgen/*.h")
file(GLOB SRC_FILES "${PROJECT_SOURCE_DIR}/include/igl/copyleft/tetgen/*.cpp")
igl_target_sources(igl_copyleft_tetgen ${INC_FILES} ${SRC_FILES})

# 4. Dependencies
include(tetgen)
igl_include(copyleft core)
target_link_libraries(igl_copyleft_tetgen ${IGL_SCOPE}
    igl::core
    igl_copyleft::core
    tetgen::tetgen
)

# 5. Unit tests
file(GLOB SRC_FILES "${PROJECT_SOURCE_DIR}/tests/include/igl/copyleft/tetgen/*.cpp")
igl_add_test(igl_copyleft_tetgen ${SRC_FILES})

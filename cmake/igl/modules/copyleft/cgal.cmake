# 1. Define module
igl_add_library(igl_copyleft_cgal)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_copyleft_cgal PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# 3. Target sources
file(GLOB INC_FILES "${PROJECT_SOURCE_DIR}/include/igl/copyleft/cgal/*.h")
file(GLOB SRC_FILES "${PROJECT_SOURCE_DIR}/include/igl/copyleft/cgal/*.cpp")
igl_target_sources(igl_copyleft_cgal ${INC_FILES} ${SRC_FILES})

# 4. Dependencies
include(cgal)
target_link_libraries(igl_copyleft_cgal ${IGL_SCOPE}
    igl::core
    igl_copyleft::core
    CGAL::CGAL
    CGAL::CGAL_Core
)

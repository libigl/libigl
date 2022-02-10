# 1. Define module
igl_add_library(igl_xml)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_xml ${IGL_SCOPE}
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# 3. Target sources
file(GLOB INC_FILES "${PROJECT_SOURCE_DIR}/include/igl/xml/*.h")
file(GLOB SRC_FILES "${PROJECT_SOURCE_DIR}/include/igl/xml/*.cpp")
igl_target_sources(igl_xml ${INC_FILES} ${SRC_FILES})

# 4. Dependencies
include(tinyxml2)
target_link_libraries(igl_xml ${IGL_SCOPE}
    igl::core
    tinyxml2::tinyxml2
)

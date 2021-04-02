# 1. Define module
igl_add_library(igl_copyleft_comiso)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_copyleft_comiso PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# 3. Target sources
file(GLOB INC_FILES "${PROJECT_SOURCE_DIR}/include/igl/copyleft/comiso/*.h")
file(GLOB SRC_FILES "${PROJECT_SOURCE_DIR}/include/igl/copyleft/comiso/*.cpp")
igl_target_sources(igl_copyleft_comiso ${INC_FILES} ${SRC_FILES})

# 4. Dependencies
include(comiso)
target_link_libraries(igl_copyleft_comiso ${IGL_SCOPE}
    igl::core
    igl_copyleft::core
    CoMISo::CoMISo
)

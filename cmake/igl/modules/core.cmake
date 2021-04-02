# 1. Define module
igl_add_library(igl_core)
set_target_properties(igl_core PROPERTIES OUTPUT_NAME igl)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_core PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# 3. Target sources
file(GLOB INC_FILES "${PROJECT_SOURCE_DIR}/include/igl/*.h")
file(GLOB SRC_FILES "${PROJECT_SOURCE_DIR}/include/igl/*.cpp")
igl_target_sources(igl_core ${INC_FILES} ${SRC_FILES})

# 4. Dependencies
include(eigen)
find_package(Threads REQUIRED)
target_link_libraries(igl_core ${IGL_SCOPE}
    Eigen3::Eigen
    Threads::Threads
)

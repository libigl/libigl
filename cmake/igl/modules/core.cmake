# 1. Define module
igl_add_library(igl_core)
if(LIBIGL_USE_STATIC_LIBRARY)
    set_target_properties(igl_core PROPERTIES OUTPUT_NAME igl)
endif()

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_core ${IGL_SCOPE}
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# 3. Target sources
file(GLOB INC_FILES "${PROJECT_SOURCE_DIR}/include/igl/*.h")
file(GLOB SRC_FILES "${PROJECT_SOURCE_DIR}/include/igl/*.cpp")
igl_target_sources(igl_core ${INC_FILES} ${SRC_FILES})

# 4. Install target & headers
igl_install(igl_core ${INC_FILES} ${SRC_FILES})

# 5. Dependencies
include(eigen)
find_package(Threads REQUIRED)
target_link_libraries(igl_core ${IGL_SCOPE}
    Eigen3::Eigen
    Threads::Threads
)

# 6. Unit tests
file(GLOB SRC_FILES "${PROJECT_SOURCE_DIR}/tests/include/igl/*.cpp")
igl_add_test(igl_core ${SRC_FILES})

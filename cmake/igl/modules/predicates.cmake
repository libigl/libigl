# 1. Define module
igl_add_library(igl_predicates)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_predicates ${IGL_SCOPE}
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# 3. Target sources
file(GLOB INC_FILES "${PROJECT_SOURCE_DIR}/include/igl/predicates/*.h")
file(GLOB SRC_FILES "${PROJECT_SOURCE_DIR}/include/igl/predicates/*.cpp")
igl_target_sources(igl_predicates ${INC_FILES} ${SRC_FILES})

# 4. Dependencies
include(predicates)
target_link_libraries(igl_predicates ${IGL_SCOPE}
    igl::core
    predicates::predicates
)

# 5. Unit tests
file(GLOB SRC_FILES "${PROJECT_SOURCE_DIR}/tests/include/igl/predicates/*.cpp")
igl_add_test(igl_predicates ${SRC_FILES})

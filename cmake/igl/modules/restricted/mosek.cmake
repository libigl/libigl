# 1. Define module
igl_add_library(igl_restricted_mosek)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_restricted_mosek ${IGL_SCOPE}
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# 3. Target sources
file(GLOB INC_FILES "${PROJECT_SOURCE_DIR}/include/igl/mosek/*.h")
file(GLOB SRC_FILES "${PROJECT_SOURCE_DIR}/include/igl/mosek/*.cpp")
igl_target_sources(igl_restricted_mosek ${INC_FILES} ${SRC_FILES})

# 4. Dependencies
find_package(MOSEK REQUIRED)
target_include_directories(igl_restricted_mosek ${IGL_SCOPE} ${MOSEK_INCLUDE_DIRS})
target_link_libraries(igl_restricted_mosek ${IGL_SCOPE} igl::core ${MOSEK_LIBRARIES})

# 5. Unit tests
file(GLOB SRC_FILES "${PROJECT_SOURCE_DIR}/tests/include/igl/mosek/*.cpp")
igl_add_test(igl_restricted_mosek ${SRC_FILES})
IF(APPLE)
  INCLUDE(${PROJECT_SOURCE_DIR}/cmake/misc/OSXFixDylibReferences.cmake)
  OSX_FIX_DYLIB_REFERENCES(test_igl_restricted_mosek "${MOSEK_LIBRARIES}")
ENDIF()

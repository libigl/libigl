# 1. Define module
igl_add_library(igl_spectra)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_spectra ${IGL_SCOPE}
    $<BUILD_INTERFACE:${libigl_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# 3. Target sources
file(GLOB INC_FILES "${libigl_SOURCE_DIR}/include/igl/spectra/*.h")
file(GLOB SRC_FILES "${libigl_SOURCE_DIR}/include/igl/spectra/*.cpp")
igl_target_sources(igl_spectra ${INC_FILES} ${SRC_FILES})

# 4. Dependencies
include(spectra)
target_link_libraries(igl_spectra ${IGL_SCOPE}
    igl::core
    spectra::spectra
)

# 5. Unit tests
file(GLOB SRC_FILES "${libigl_SOURCE_DIR}/tests/include/igl/spectra/*.cpp")
igl_add_test(igl_spectra ${SRC_FILES})

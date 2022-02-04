# 1. Define module
igl_add_library(igl_restricted_matlab)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_restricted_matlab ${IGL_SCOPE}
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# 3. Target sources
file(GLOB INC_FILES "${PROJECT_SOURCE_DIR}/include/igl/matlab/*.h")
file(GLOB SRC_FILES "${PROJECT_SOURCE_DIR}/include/igl/matlab/*.cpp")
igl_target_sources(igl_restricted_matlab ${INC_FILES} ${SRC_FILES})

# 4. Dependencies
# MATLAB_ADDITIONAL_VERSIONS should have already been set in libigl/CMakeLists.txt
find_package(Matlab REQUIRED COMPONENTS MEX_COMPILER MX_LIBRARY ENG_LIBRARY MAT_LIBRARY)
# to-do detect if imported targets are available (cmake 3.22's findmatlab has
# them)  and use those instead of hard coding these.
target_link_libraries(igl_restricted_matlab ${IGL_SCOPE} igl::core ${Matlab_LIBRARIES})
target_include_directories(igl_restricted_matlab ${IGL_SCOPE} ${Matlab_INCLUDE_DIRS})

# 5. Unit tests
# to-do write some matlab tests

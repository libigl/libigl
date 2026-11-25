# 1. Define module
igl_add_library(igl_core)
if(LIBIGL_USE_STATIC_LIBRARY)
    set_target_properties(igl_core PROPERTIES OUTPUT_NAME igl)
endif()

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_core ${IGL_SCOPE}
    $<BUILD_INTERFACE:${libigl_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# 3. Target sources
file(GLOB INC_FILES "${libigl_SOURCE_DIR}/include/igl/*.h")
file(GLOB SRC_FILES "${libigl_SOURCE_DIR}/include/igl/*.cpp")
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

# 5b. igl::parallel_for backend. parallel_for.h is header-only, so the selected
# backend's definition (and any compile flags / link) must propagate to every
# consumer — hence ${IGL_SCOPE} (PUBLIC for static, INTERFACE for header-only).
# POOL is the dependency-free default; a requested-but-missing OPENMP/TBB warns
# and falls back to POOL so the build never breaks on a missing optional dep.
if(NOT DEFINED LIBIGL_PARALLEL_FOR_BACKEND OR LIBIGL_PARALLEL_FOR_BACKEND STREQUAL "")
    set(LIBIGL_PARALLEL_FOR_BACKEND "POOL")
endif()
if(LIBIGL_PARALLEL_FOR_BACKEND STREQUAL "POOL")
    # internal std::thread pool: nothing to add
elseif(LIBIGL_PARALLEL_FOR_BACKEND STREQUAL "SERIAL")
    target_compile_definitions(igl_core ${IGL_SCOPE} IGL_PARALLEL_FOR_FORCE_SERIAL)
elseif(LIBIGL_PARALLEL_FOR_BACKEND STREQUAL "OPENMP")
    find_package(OpenMP QUIET COMPONENTS CXX)
    if(OpenMP_CXX_FOUND)
        message(STATUS "igl::parallel_for backend: OpenMP")
        target_compile_definitions(igl_core ${IGL_SCOPE} IGL_PARALLEL_FOR_OPENMP)
        target_link_libraries(igl_core ${IGL_SCOPE} OpenMP::OpenMP_CXX)
    else()
        message(WARNING "LIBIGL_PARALLEL_FOR_BACKEND=OPENMP but OpenMP was not found; falling back to the internal thread pool (POOL).")
    endif()
elseif(LIBIGL_PARALLEL_FOR_BACKEND STREQUAL "TBB")
    find_package(TBB QUIET)
    if(TARGET TBB::tbb)
        message(STATUS "igl::parallel_for backend: TBB")
        target_compile_definitions(igl_core ${IGL_SCOPE} IGL_PARALLEL_FOR_TBB)
        target_link_libraries(igl_core ${IGL_SCOPE} TBB::tbb)
    else()
        message(WARNING "LIBIGL_PARALLEL_FOR_BACKEND=TBB but TBB was not found; falling back to the internal thread pool (POOL).")
    endif()
else()
    message(FATAL_ERROR "Unknown LIBIGL_PARALLEL_FOR_BACKEND '${LIBIGL_PARALLEL_FOR_BACKEND}' (expected POOL, OPENMP, TBB, or SERIAL).")
endif()

# 6. Unit tests
file(GLOB SRC_FILES "${libigl_SOURCE_DIR}/tests/include/igl/*.cpp")
igl_add_test(igl_core ${SRC_FILES})

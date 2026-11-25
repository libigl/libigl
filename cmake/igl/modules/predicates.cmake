# 1. Define module
igl_add_library(igl_predicates)

# 2. Include headers
include(GNUInstallDirs)
target_include_directories(igl_predicates ${IGL_SCOPE}
    $<BUILD_INTERFACE:${libigl_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# 3. Target sources
file(GLOB INC_FILES "${libigl_SOURCE_DIR}/include/igl/predicates/*.h")
file(GLOB SRC_FILES "${libigl_SOURCE_DIR}/include/igl/predicates/*.cpp")
igl_target_sources(igl_predicates ${INC_FILES} ${SRC_FILES})

# 4. Dependencies
include(predicates)
target_link_libraries(igl_predicates ${IGL_SCOPE}
    igl::core
    predicates::predicates
)

# Preserve the rounded operations used by the exact sign filters, including
# when their implementations are compiled in a header-only consumer.
if(MSVC)
    target_compile_options(igl_predicates ${IGL_SCOPE} /fp:strict)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(igl_predicates ${IGL_SCOPE} -fno-fast-math -ffp-contract=off)
endif()

# 5. Unit tests
file(GLOB SRC_FILES "${libigl_SOURCE_DIR}/tests/include/igl/predicates/*.cpp")
igl_add_test(igl_predicates ${SRC_FILES})
if(LIBIGL_BUILD_TESTS)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        add_executable(test_predicates_allocations "${libigl_SOURCE_DIR}/tests/predicates_allocations.cpp")
        target_link_libraries(test_predicates_allocations PRIVATE igl::predicates)
        target_link_options(test_predicates_allocations PRIVATE
            -Wl,--wrap=predicates_expansion_sum -Wl,--wrap=nextafter)
        add_test(NAME predicates_allocations COMMAND test_predicates_allocations)
        set_tests_properties(predicates_allocations PROPERTIES TIMEOUT 30)
    endif()
    get_target_property(predicates_source_dir predicates SOURCE_DIR)
    add_test(NAME predicates_patch
        COMMAND ${CMAKE_COMMAND}
            "-DPREDICATES_SOURCE_DIR=${predicates_source_dir}"
            "-DPATCH_DIR=${libigl_SOURCE_DIR}/cmake/recipes/external"
            "-DTEST_DIR=${CMAKE_CURRENT_BINARY_DIR}/predicates-patch-test"
            -P "${libigl_SOURCE_DIR}/tests/cmake/predicates_patch.cmake")
    set_tests_properties(predicates_patch PROPERTIES TIMEOUT 60)
endif()

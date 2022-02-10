function(igl_add_test module_name)
    if(NOT LIBIGL_BUILD_TESTS)
        return()
    endif()

    if(NOT TARGET ${module_name})
        message(FATAL_ERROR "'${module_name}' is not a CMake target")
    endif()

    # Check if category is `copyleft` or `restricted`
    if(${module_name} MATCHES "^igl_copyleft")
        set(suffix "_copyleft")
    elseif(${module_name} MATCHES "^igl_restricted")
        set(suffix "_restricted")
    else()
        set(suffix "")
    endif()

    # Create test executable
    add_executable(test_${module_name}
        ${PROJECT_SOURCE_DIR}/tests/main.cpp
        ${PROJECT_SOURCE_DIR}/tests/test_common.h
        ${ARGN}
    )

    # Include headers
    target_include_directories(test_${module_name} PUBLIC ${PROJECT_SOURCE_DIR}/tests)

    # Compile definitions
    target_compile_definitions(test_${module_name} PUBLIC CATCH_CONFIG_ENABLE_BENCHMARKING)

    # Dependencies
    include(catch2)
    include(libigl_tests_data)
    target_link_libraries(test_${module_name} PUBLIC
        ${module_name}
        igl::tests_data
        Catch2::Catch2
    )

    # IDE Folder
    set_target_properties(test_${module_name} PROPERTIES FOLDER Libigl_Tests)

    # Output directory
    set_target_properties(test_${module_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests")

    # Register tests
    FetchContent_GetProperties(catch2)
    include("${catch2_SOURCE_DIR}/contrib/Catch.cmake")
    catch_discover_tests(test_${module_name})
endfunction()

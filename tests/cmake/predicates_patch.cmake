# Exercise dependency patching with both Unix and Windows line endings.
cmake_minimum_required(VERSION 3.16)
find_package(Git REQUIRED)
foreach(required PREDICATES_SOURCE_DIR PATCH_DIR TEST_DIR)
    if(NOT DEFINED ${required} OR "${${required}}" STREQUAL "")
        message(FATAL_ERROR "Set ${required}")
    endif()
endforeach()
function(run)
    execute_process(COMMAND ${ARGV} RESULT_VARIABLE result
        OUTPUT_VARIABLE output ERROR_VARIABLE error TIMEOUT 10)
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "Command failed (${result}): ${ARGV}\n${output}\n${error}")
    endif()
endfunction()
function(check_newline_style path style)
    # Inspect the actual bytes, including any doubled carriage returns.
    file(READ "${path}" bytes HEX)
    if(style STREQUAL CRLF)
        if(NOT bytes MATCHES "0d0a")
            message(FATAL_ERROR "${path} has no CRLF line endings")
        endif()
        string(REPLACE "0d0a" "" bytes "${bytes}")
        if(bytes MATCHES "0a")
            message(FATAL_ERROR "${path} contains a bare LF")
        endif()
    endif()
    if(bytes MATCHES "0d")
        message(FATAL_ERROR "${path} contains an unexpected carriage return")
    endif()
endfunction()
file(MAKE_DIRECTORY "${TEST_DIR}")
# file(WRITE) uses text mode and converts LF back to CRLF on Windows.
# Specify the on-disk line endings for every source and patch fixture instead.
configure_file("${PATCH_DIR}/predicates-expansion.patch"
    "${TEST_DIR}/normalized.patch" @ONLY NEWLINE_STYLE UNIX)
check_newline_style("${TEST_DIR}/normalized.patch" LF)
set(files CMakeLists.txt predicates.c predicates.h)
foreach(source_eol LF CRLF)
    foreach(patch_eol LF CRLF)
        set(fixture "${TEST_DIR}/${source_eol}-${patch_eol}")
        file(MAKE_DIRECTORY "${fixture}/source" "${fixture}/patch")
        # A separate Git root ensures git apply cannot silently skip paths
        # because this fixture is nested inside the libigl checkout.
        run("${GIT_EXECUTABLE}" init --quiet "${fixture}/source")
        foreach(name IN LISTS files)
            configure_file("${PREDICATES_SOURCE_DIR}/${name}"
                "${fixture}/source/${name}" @ONLY NEWLINE_STYLE UNIX)
        endforeach()
        run("${GIT_EXECUTABLE}" -C "${fixture}/source" apply --reverse
            --ignore-space-change "${TEST_DIR}/normalized.patch")
        # git apply may follow core.autocrlf when restoring the source, so
        # enforce the requested fixture style after reversing the patch too.
        foreach(name IN LISTS files)
            configure_file("${fixture}/source/${name}"
                "${fixture}/${name}" @ONLY NEWLINE_STYLE "${source_eol}")
            file(RENAME "${fixture}/${name}" "${fixture}/source/${name}")
            check_newline_style("${fixture}/source/${name}" "${source_eol}")
        endforeach()
        configure_file("${TEST_DIR}/normalized.patch"
            "${fixture}/patch/predicates-expansion.patch" @ONLY NEWLINE_STYLE "${patch_eol}")
        check_newline_style("${fixture}/patch/predicates-expansion.patch" "${patch_eol}")
        configure_file("${PATCH_DIR}/predicates-patch.cmake"
            "${fixture}/patch/predicates-patch.cmake" COPYONLY)
        foreach(pass RANGE 1 2)
            run("${CMAKE_COMMAND}" "-DSOURCE_DIR=${fixture}/source"
                -P "${fixture}/patch/predicates-patch.cmake")
            foreach(name IN LISTS files)
                file(READ "${PREDICATES_SOURCE_DIR}/${name}" expected)
                file(READ "${fixture}/source/${name}" actual)
                string(REPLACE "\r\n" "\n" expected "${expected}")
                string(REPLACE "\r\n" "\n" actual "${actual}")
                if(NOT actual STREQUAL expected)
                    message(FATAL_ERROR "${source_eol}/${patch_eol} pass ${pass} changed ${name} incorrectly")
                endif()
            endforeach()
        endforeach()
        # Do not accept an unexpected source tree as already patched.
        file(WRITE "${fixture}/source/predicates.h" "unexpected dependency version\n")
        execute_process(COMMAND "${CMAKE_COMMAND}" "-DSOURCE_DIR=${fixture}/source"
            -P "${fixture}/patch/predicates-patch.cmake" RESULT_VARIABLE rejected
            OUTPUT_QUIET ERROR_VARIABLE rejection TIMEOUT 10)
        if(rejected EQUAL 0 OR NOT rejection MATCHES "Predicates source does not match")
            message(FATAL_ERROR "${source_eol}/${patch_eol} did not reject unexpected source")
        endif()
        message(STATUS "${source_eol}/${patch_eol}: apply, repeat and source rejection passed")
    endforeach()
endforeach()
file(REMOVE_RECURSE "${TEST_DIR}")

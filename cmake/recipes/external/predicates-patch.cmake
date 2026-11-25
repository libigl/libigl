find_package(Git REQUIRED)
# Windows checkouts may convert both the patch and its target files to CRLF.
# Git cannot parse CRLF patch files, even when it ignores context whitespace.
set(patch "${SOURCE_DIR}/libigl-predicates-expansion.patch")
configure_file("${CMAKE_CURRENT_LIST_DIR}/predicates-expansion.patch" "${patch}"
    @ONLY NEWLINE_STYLE UNIX)
execute_process(COMMAND "${GIT_EXECUTABLE}" apply --ignore-space-change --check "${patch}"
    WORKING_DIRECTORY "${SOURCE_DIR}" RESULT_VARIABLE can_apply
    OUTPUT_QUIET ERROR_VARIABLE forward_error)
if(can_apply EQUAL 0)
    execute_process(COMMAND "${GIT_EXECUTABLE}" apply --ignore-space-change "${patch}"
        WORKING_DIRECTORY "${SOURCE_DIR}" RESULT_VARIABLE applied
        OUTPUT_QUIET ERROR_VARIABLE apply_error)
    if(NOT applied EQUAL 0)
        set(patch_error "Could not expose predicates expansion arithmetic: ${apply_error}")
    endif()
else()
    # FetchContent may repeat its patch step on a subsequent configure.
    execute_process(COMMAND "${GIT_EXECUTABLE}" apply --ignore-space-change --reverse --check "${patch}"
        WORKING_DIRECTORY "${SOURCE_DIR}" RESULT_VARIABLE already_applied
        OUTPUT_QUIET ERROR_VARIABLE reverse_error)
    if(NOT already_applied EQUAL 0)
        set(patch_error "Predicates source does not match the expansion arithmetic patch.\nForward check: ${forward_error}\nReverse check: ${reverse_error}")
    endif()
endif()
file(REMOVE "${patch}")
if(DEFINED patch_error)
    message(FATAL_ERROR "${patch_error}")
endif()

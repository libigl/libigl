find_package(Git REQUIRED)
set(patch "${CMAKE_CURRENT_LIST_DIR}/predicates-expansion.patch")
execute_process(COMMAND "${GIT_EXECUTABLE}" apply --check "${patch}"
    WORKING_DIRECTORY "${SOURCE_DIR}" RESULT_VARIABLE can_apply
    OUTPUT_QUIET ERROR_QUIET)
if(can_apply EQUAL 0)
    execute_process(COMMAND "${GIT_EXECUTABLE}" apply "${patch}"
        WORKING_DIRECTORY "${SOURCE_DIR}" RESULT_VARIABLE applied)
    if(NOT applied EQUAL 0)
        message(FATAL_ERROR "Could not expose predicates expansion arithmetic")
    endif()
else()
    # FetchContent may repeat its patch step on a subsequent configure.
    execute_process(COMMAND "${GIT_EXECUTABLE}" apply --reverse --check "${patch}"
        WORKING_DIRECTORY "${SOURCE_DIR}" RESULT_VARIABLE already_applied
        OUTPUT_QUIET ERROR_QUIET)
    if(NOT already_applied EQUAL 0)
        message(FATAL_ERROR "Predicates source does not match the expansion arithmetic patch")
    endif()
endif()

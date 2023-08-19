# Transitively list all link libraries of a target (recursive call)
function(igl_get_dependencies_recursive OUTPUT_VARIABLE TARGET)
    get_target_property(_aliased ${TARGET} ALIASED_TARGET)
    if(_aliased)
        set(TARGET ${_aliased})
    endif()

    get_target_property(_IMPORTED ${TARGET} IMPORTED)
    get_target_property(_TYPE ${TARGET} TYPE)
    if(_IMPORTED OR (${_TYPE} STREQUAL "INTERFACE_LIBRARY"))
        get_target_property(TARGET_DEPENDENCIES ${TARGET} INTERFACE_LINK_LIBRARIES)
    else()
        get_target_property(TARGET_DEPENDENCIES ${TARGET} LINK_LIBRARIES)
    endif()

    # MKL-specific list of runtime dependencies
    get_property(RUNTIME_DEPENDENCIES TARGET ${TARGET} PROPERTY mkl_RUNTIME_DEPENDENCIES)
    if(RUNTIME_DEPENDENCIES)
        list(APPEND TARGET_DEPENDENCIES ${RUNTIME_DEPENDENCIES})
    endif()

    set(VISITED_TARGETS ${${OUTPUT_VARIABLE}})
    foreach(DEPENDENCY IN ITEMS ${TARGET_DEPENDENCIES})
        if(TARGET ${DEPENDENCY})
            get_target_property(_aliased ${DEPENDENCY} ALIASED_TARGET)
            if(_aliased)
                set(DEPENDENCY ${_aliased})
            endif()

            if(NOT (DEPENDENCY IN_LIST VISITED_TARGETS))
                list(APPEND VISITED_TARGETS ${DEPENDENCY})
                igl_get_dependencies_recursive(VISITED_TARGETS ${DEPENDENCY})
            endif()
        endif()
    endforeach()
    set(${OUTPUT_VARIABLE} ${VISITED_TARGETS} PARENT_SCOPE)
endfunction()

# Transitively list all link libraries of a target
function(igl_get_dependencies OUTPUT_VARIABLE TARGET)
    set(DISCOVERED_TARGETS "")
    igl_get_dependencies_recursive(DISCOVERED_TARGETS ${TARGET})
    set(${OUTPUT_VARIABLE} ${DISCOVERED_TARGETS} PARENT_SCOPE)
endfunction()

# Copy .dll dependencies to a target executable's folder. This function must be called *after* all the CMake
# dependencies of the executable target have been defined, otherwise some .dlls might not be copied to the target
# folder.
function(igl_copy_dll target)
    if(NOT WIN32)
        return()
    endif()

    # Sanity checks
    get_target_property(TYPE ${target} TYPE)
    if(NOT ${TYPE} STREQUAL "EXECUTABLE")
        message(FATAL_ERROR "igl_copy_dll() was called on a non-executable target: ${target}")
    endif()

    # Create a custom command to do the actual copy. This needs to be executed before Catch2's POST_BUILD command,
    # so we define this as a PRE_LINK command for the executable target.
    add_custom_command(
        TARGET ${target}
        PRE_LINK
        COMMAND ${CMAKE_COMMAND} -E touch "${CMAKE_BINARY_DIR}/runtime_deps/copy_dll_${target}_$<CONFIG>.cmake"
        COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/runtime_deps/copy_dll_${target}_$<CONFIG>.cmake"
        COMMENT "Copying dlls for target ${target}"
    )

    # Retrieve all target dependencies
    igl_get_dependencies(TARGET_DEPENDENCIES ${target})

    # Iterate over dependencies, and create a copy rule for each .dll that we find
    set(COPY_SCRIPT_CONTENT "")
    foreach(DEPENDENCY IN ITEMS ${TARGET_DEPENDENCIES})
        get_target_property(TYPE ${DEPENDENCY} TYPE)
        if(NOT (${TYPE} STREQUAL "SHARED_LIBRARY" OR ${TYPE} STREQUAL "MODULE_LIBRARY"))
            continue()
        endif()

        # Instruction to copy target file if it exists
        string(APPEND COPY_SCRIPT_CONTENT
            "if(EXISTS \"$<TARGET_FILE:${DEPENDENCY}>\")\n    "
                "execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different "
                "\"$<TARGET_FILE:${DEPENDENCY}>\" "
                "\"$<TARGET_FILE_DIR:${target}>/$<TARGET_FILE_NAME:${DEPENDENCY}>\")\n"
            "endif()\n"
        )
    endforeach()

    # Finally generate one script for each configuration supported by this generator
    message(STATUS "Populating copy rules for target: ${target}")
    file(GENERATE
        OUTPUT ${CMAKE_BINARY_DIR}/runtime_deps/copy_dll_${target}_$<CONFIG>.cmake
        CONTENT "${COPY_SCRIPT_CONTENT}"
    )
endfunction()

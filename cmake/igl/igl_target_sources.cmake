function(igl_target_sources module_name)
    if(${CMAKE_VERSION} VERSION_LESS 3.19 AND NOT LIBIGL_USE_STATIC_LIBRARY)
        # Old approach: defines a dummy custom target for the IDE
        add_custom_target(${module_name}_ SOURCES ${ARGN})
    else()
        target_sources(${module_name} PRIVATE ${ARGN})
    endif()
endfunction()

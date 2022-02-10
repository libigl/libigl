# Helper functions to include libigl modules
function(_igl_include_full prefix name force)
    string(TOUPPER "${prefix}" prefix_uc)
    string(TOUPPER "${name}" name_uc)
    if(prefix_uc)
        string(PREPEND prefix_uc _)
    endif()
    string(TOLOWER "${prefix_uc}" prefix_lc)

    if(TARGET igl${prefix_lc}::${name})
        # Target already exists, skip
        return()
    endif()

    if(${force} AND NOT ${name} STREQUAL "core")
        message(STATUS "Forcing include of libigl module: ${name}")
    endif()

    # Include igl target definition
    if(LIBIGL${prefix_uc}_WITH_${name_uc} OR ${force})
        include(${PROJECT_SOURCE_DIR}/cmake/igl/modules/${prefix}/${name}.cmake)
    endif()
endfunction()

# Include module only if CMake option is provided
function(igl_include_optional name)
    if(ARGC GREATER_EQUAL 2)
        # Two args given: prefix + name
        _igl_include_full(${name} ${ARGV1} FALSE)
    else()
        # Only one arg given: prefix (which contains module name)
        _igl_include_full("" ${name} FALSE)
    endif()
endfunction()

# Include module only unconditionally (e.g. if module is a dependency of another module)
function(igl_include name)
    if(ARGC GREATER_EQUAL 2)
        # Two args given: prefix + name
        _igl_include_full(${name} ${ARGV1} TRUE)
    else()
        # Only one arg given: prefix (which contains module name)
        _igl_include_full("" ${name} TRUE)
    endif()
endfunction()

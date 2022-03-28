# Try to find the GNU Multiple Precision Arithmetic Library (GMP)
# See http://gmplib.org/

if(${CMAKE_VERSION} VERSION_LESS "3.18.0")
    set(REQUIRED_FLAG "")
else()
    set(REQUIRED_FLAG REQUIRED)
endif()

# On Windows, we must use the pre-compiled versions downloaded with libigl
if(WIN32)
    set(NO_DEFAULT_FLAG NO_DEFAULT_PATH)
else()
    set(NO_DEFAULT_FLAG "")
endif()

find_path(GMP_INCLUDES
    NAMES
        gmp.h
    PATHS
        ENV GMP_DIR
        ${INCLUDE_INSTALL_DIR}
    PATH_SUFFIXES
        include
    ${REQUIRED_FLAG}
    ${NO_DEFAULT_FLAG}
)

find_library(GMP_LIBRARIES
    NAMES
        gmp
        libgmp-10
    PATHS
        ENV GMP_DIR
        ${LIB_INSTALL_DIR}
    PATH_SUFFIXES
        lib
    ${REQUIRED_FLAG}
    ${NO_DEFAULT_FLAG}
)

set(GMP_EXTRA_VARS "")
if(WIN32)
    # Find dll file and set IMPORTED_LOCATION to the .dll file
    find_file(GMP_RUNTIME_LIB
        NAMES
            gmp.dll
            libgmp-10.dll
        PATHS
            ENV GMP_DIR
            ${LIB_INSTALL_DIR}
        PATH_SUFFIXES
            lib
        ${REQUIRED_FLAG}
        ${NO_DEFAULT_FLAG}
    )
    list(APPEND GMP_EXTRA_VARS GMP_RUNTIME_LIB)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP
    REQUIRED_VARS
        GMP_INCLUDES
        GMP_LIBRARIES
        ${GMP_EXTRA_VARS}
    REASON_FAILURE_MESSAGE
        "GMP is not installed on your system. Either install GMP using your preferred package manager, or disable libigl modules that depend on GMP, such as CGAL. See LibiglOptions.cmake.sample for configuration options. Do not forget to delete your <build>/CMakeCache.txt for the changes to take effect."
)
mark_as_advanced(GMP_INCLUDES GMP_LIBRARIES)

if(GMP_INCLUDES AND GMP_LIBRARIES AND NOT TARGET gmp::gmp)
    if(GMP_RUNTIME_LIB)
        add_library(gmp::gmp SHARED IMPORTED)
    else()
        add_library(gmp::gmp UNKNOWN IMPORTED)
    endif()

    # Set public header location and link language
    set_target_properties(gmp::gmp PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "C"
        INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDES}"
    )

    # Set lib location. On Windows we specify both the .lib and the .dll paths
    if(GMP_RUNTIME_LIB)
        set_target_properties(gmp::gmp PROPERTIES
            IMPORTED_IMPLIB "${GMP_LIBRARIES}"
            IMPORTED_LOCATION "${GMP_RUNTIME_LIB}"
        )
    else()
        set_target_properties(gmp::gmp PROPERTIES
            IMPORTED_LOCATION "${GMP_LIBRARIES}"
        )
    endif()
endif()

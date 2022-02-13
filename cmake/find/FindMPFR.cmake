# Try to find the MPFR library
# See http://www.mpfr.org/

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

find_path(MPFR_INCLUDES
    NAMES
        mpfr.h
    PATHS
        ENV MPFR_DIR
        ${INCLUDE_INSTALL_DIR}
    PATH_SUFFIXES
        include
    ${REQUIRED_FLAG}
    ${NO_DEFAULT_FLAG}
)

find_library(MPFR_LIBRARIES
    NAMES
        mpfr
        libmpfr-4
    PATHS
        ENV MPFR_DIR
        ${LIB_INSTALL_DIR}
    PATH_SUFFIXES
        lib
    ${REQUIRED_FLAG}
    ${NO_DEFAULT_FLAG}
)

set(MPFR_EXTRA_VARS "")
if(WIN32)
    # Find dll file and set IMPORTED_LOCATION to the .dll file
    find_file(MPFR_RUNTIME_LIB
        NAMES
            mpfr.dll
            libmpfr-4.dll
        PATHS
            ENV MPFR_DIR
            ${LIB_INSTALL_DIR}
        PATH_SUFFIXES
            lib
        ${REQUIRED_FLAG}
        ${NO_DEFAULT_FLAG}
    )
    list(APPEND MPFR_EXTRA_VARS MPFR_RUNTIME_LIB)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPFR
    REQUIRED_VARS
        MPFR_INCLUDES
        MPFR_LIBRARIES
        ${MPFR_EXTRA_VARS}
    REASON_FAILURE_MESSAGE
        "MPFR is not installed on your system. Either install MPFR using your preferred package manager, or disable libigl modules that depend on MPFR, such as CGAL. See LibiglOptions.cmake.sample for configuration options. Do not forget to delete your <build>/CMakeCache.txt for the changes to take effect."
)
mark_as_advanced(MPFR_INCLUDES MPFR_LIBRARIES)

if(MPFR_INCLUDES AND MPFR_LIBRARIES AND NOT TARGET mpfr::mpfr)
    if(MPFR_RUNTIME_LIB)
        add_library(mpfr::mpfr SHARED IMPORTED)
    else()
        add_library(mpfr::mpfr UNKNOWN IMPORTED)
    endif()

    # Set public header location and link language
    set_target_properties(mpfr::mpfr PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "C"
        INTERFACE_INCLUDE_DIRECTORIES "${MPFR_INCLUDES}"
    )

    # Set lib location. On Windows we specify both the .lib and the .dll paths
    if(MPFR_RUNTIME_LIB)
        set_target_properties(mpfr::mpfr PROPERTIES
            IMPORTED_IMPLIB "${MPFR_LIBRARIES}"
            IMPORTED_LOCATION "${MPFR_RUNTIME_LIB}"
        )
    else()
        set_target_properties(mpfr::mpfr PROPERTIES
            IMPORTED_LOCATION "${MPFR_LIBRARIES}"
        )
    endif()
endif()

# Try to find the MPFR library
# See http://www.mpfr.org/

find_path(MPFR_INCLUDES
    NAMES
        mpfr.h
    PATHS
        ENV MPFR_DIR
        ${INCLUDE_INSTALL_DIR}
    PATH_SUFFIXES
        include
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
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPFR
    REQUIRED_VARS
        MPFR_INCLUDES
        MPFR_LIBRARIES
    REASON_FAILURE_MESSAGE
        "MPFR is not installed on your system. Either install MPFR using your preferred package manager, or disable libigl modules that depend on MPFR, such as CGAL. See LibiglOptions.cmake.sample for configuration options. Do not forget to delete your <build>/CMakeCache.txt for the changes to take effect."
)
mark_as_advanced(MPFR_INCLUDES MPFR_LIBRARIES)

if(MPFR_INCLUDES AND MPFR_LIBRARIES AND NOT TARGET mpfr::mpfr)
    add_library(mpfr::mpfr UNKNOWN IMPORTED)
    set_target_properties(mpfr::mpfr PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "C"
        IMPORTED_LOCATION "${MPFR_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${MPFR_INCLUDES}"
    )
endif()

# Try to find the MPFR library
# See http://www.mpfr.org/

# Note: We use MPFR_INCLUDE_DIR and MPFR_LIBRARIES to match with CGAL's FindMPFR.cmake module.

find_path(MPFR_INCLUDE_DIR
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
        MPFR_INCLUDE_DIR
        MPFR_LIBRARIES
    REASON_FAILURE_MESSAGE
        "MPFR is not installed on your system. Either install MPFR using your preferred package manager, or disable libigl modules that depend on MPFR, such as CGAL. See LibiglOptions.cmake.sample for configuration options. Do not forget to delete your <build>/CMakeCache.txt for the changes to take effect."
)
mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIBRARIES)

if(MPFR_INCLUDE_DIR AND MPFR_LIBRARIES AND NOT TARGET mpfr::mpfr)
    add_library(mpfr::mpfr UNKNOWN IMPORTED)
    set_target_properties(mpfr::mpfr PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "C"
        IMPORTED_LOCATION "${MPFR_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${MPFR_INCLUDE_DIR}"
    )
endif()

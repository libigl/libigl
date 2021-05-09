# Try to find the GNU Multiple Precision Arithmetic Library (GMP)
# See http://gmplib.org/

# Note: We use GMP_INCLUDE_DIR and GMP_LIBRARIES to match with CGAL's FindGMP.cmake module.

find_path(GMP_INCLUDE_DIR
    NAMES
        gmp.h
    PATHS
        ENV GMP_DIR
        ${INCLUDE_INSTALL_DIR}
    PATH_SUFFIXES
        include
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
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP
    REQUIRED_VARS
        GMP_INCLUDE_DIR
        GMP_LIBRARIES
    REASON_FAILURE_MESSAGE
        "GMP is not installed on your system. Either install GMP using your preferred package manager, or disable libigl modules that depend on GMP, such as CORK and CGAL. See LibiglOptions.cmake.sample for configuration options. Do not forget to delete your <build>/CMakeCache.txt for the changes to take effect."
)
mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARIES)

if(GMP_INCLUDE_DIR AND GMP_LIBRARIES AND NOT TARGET gmp::gmp)
    add_library(gmp::gmp UNKNOWN IMPORTED)
    set_target_properties(gmp::gmp PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "C"
        IMPORTED_LOCATION "${GMP_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}"
    )
endif()

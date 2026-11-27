# Finds the Eigen installed alongside libigl (see CMakeLists.txt).
if(NOT TARGET Eigen3::Eigen)
    find_file(EIGEN_TARGETS_FILE EigenTargets.cmake
        HINTS "${libigl_DIR}/../eigen"
        PATH_SUFFIXES lib/cmake/eigen
    )
    if(EIGEN_TARGETS_FILE)
        include("${EIGEN_TARGETS_FILE}")
    endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Eigen3 REQUIRED_VARS EIGEN_TARGETS_FILE)

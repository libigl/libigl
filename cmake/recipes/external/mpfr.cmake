if(TARGET mpfr::mpfr)
    return()
endif()

# Download precompiled .dll on Windows
include(gmp_mpfr)

message(STATUS "Third-party: creating target 'mpfr::mpfr'")

# Find_package will look for our downloaded lib on Windows, and system-wide on Linux/macOS
find_package(MPFR REQUIRED)

if(NOT TARGET mpfr::mpfr)
    message(FATAL_ERROR "Creation of target 'mpfr::mpfr' failed")
endif()

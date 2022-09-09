if(TARGET gmp::gmp)
    return()
endif()

# Download precompiled .dll on Windows
include(gmp_mpfr)

message(STATUS "Third-party: creating target 'gmp::gmp'")

# Find_package will look for our downloaded lib on Windows, and system-wide on Linux/macOS
find_package(GMP REQUIRED)

if(NOT TARGET gmp::gmp)
    message(FATAL_ERROR "Creation of target 'gmp::gmp' failed")
endif()

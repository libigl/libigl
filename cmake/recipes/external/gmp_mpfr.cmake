if(WIN32)
    message(STATUS "Third-party: downloading gmp + mpfr")

    include(FetchContent)
    FetchContent_Declare(
        gmp_mpfr
        URL https://github.com/CGAL/cgal/releases/download/v5.2.1/CGAL-5.2.1-win64-auxiliary-libraries-gmp-mpfr.zip
        URL_MD5 247f4dca741c6b9a9be76286414070fa
    )
    FetchContent_MakeAvailable(gmp_mpfr)

    # For CGAL and Cork
    set(ENV{GMP_DIR} "${gmp_mpfr_SOURCE_DIR}/gmp")
    set(ENV{MPFR_DIR} "${gmp_mpfr_SOURCE_DIR}/gmp")
else()
    # On Linux/macOS, gmp+mpfr should be installed system-wide
endif()

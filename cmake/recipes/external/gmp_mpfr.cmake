if(WIN32)
    message(STATUS "Third-party: downloading gmp + mpfr")

    include(FetchContent)

    # CGAL 5+ ships with a single .zip combining GMP + MPFR's precompiled dlls.
    # For now we still download them separately.

    FetchContent_Declare(
        gmp
        URL     https://cgal.geometryfactory.com/CGAL/precompiled_libs/auxiliary/x64/GMP/5.0.1/gmp-all-CGAL-3.9.zip
        URL_MD5 508c1292319c832609329116a8234c9f
    )
    FetchContent_MakeAvailable(gmp)

    FetchContent_Declare(
        mpfr
        URL https://cgal.geometryfactory.com/CGAL/precompiled_libs/auxiliary/x64/MPFR/3.0.0/mpfr-all-CGAL-3.9.zip
        URL_MD5 48840454eef0ff18730050c05028734b
    )
    FetchContent_MakeAvailable(mpfr)

    # FetchContent_Declare(
    #     gmp_mpfr
    #     URL https://github.com/CGAL/cgal/releases/download/v5.2.1/CGAL-5.2.1-win64-auxiliary-libraries-gmp-mpfr.zip
    #     URL_MD5 247f4dca741c6b9a9be76286414070fa
    # )

    # For CGAL
    set(ENV{GMP_DIR} "${gmp_SOURCE_DIR}")
    set(ENV{MPFR_DIR} "${mpfr_SOURCE_DIR}")
else()
    # On Linux/macOS, gmp+mpfr should be installed system-wide
endif()

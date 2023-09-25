if(TARGET CGAL::CGAL)
    return()
endif()

message(STATUS "Third-party: creating target 'CGAL::CGAL'")

include(FetchContent)
FetchContent_Declare(
    cgal
    URL https://github.com/CGAL/cgal/releases/download/v5.6/CGAL-5.6-library.tar.xz
    URL_MD5 793da2d1597f3a5c0e3524f73a0b4039
)
FetchContent_GetProperties(cgal)
if(cgal_POPULATED)
    return()
endif()
FetchContent_Populate(cgal)

function(cgal_import_target)
    macro(ignore_package NAME VERSION_NUM)
        include(CMakePackageConfigHelpers)
        file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${NAME}/${NAME}Config.cmake "")
        write_basic_package_version_file(
            ${CMAKE_CURRENT_BINARY_DIR}/${NAME}/${NAME}ConfigVersion.cmake
            VERSION ${VERSION_NUM}
            COMPATIBILITY AnyNewerVersion
        )
        set(${NAME}_DIR ${CMAKE_CURRENT_BINARY_DIR}/${NAME} CACHE PATH "")
        set(${NAME}_ROOT ${CMAKE_CURRENT_BINARY_DIR}/${NAME} CACHE PATH "")
    endmacro()

    include(gmp)
    include(mpfr)
    include(boost)

    ignore_package(GMP 5.0.1)
    set(GMP_INCLUDE_DIR ${gmp_INCLUDE_DIR})
    set(GMP_LIBRARIES gmp::gmp)
    set(GMPXX_INCLUDE_DIR ${GMP_INCLUDE_DIR})
    set(GMPXX_LIBRARIES ${GMP_LIBRARIES})

    ignore_package(MPFR 3.0.0)
    set(MPFR_INCLUDE_DIR "")
    set(MPFR_LIBRARIES mpfr::mpfr)

    ignore_package(Boost 1.71.0)
    set(Boost_INCLUDE_DIRS "")
    set(Boost_LIBRARIES Boost::thread Boost::system)

    # Prefer Config mode before Module mode to prevent CGAL from loading its own FindXXX.cmake
    set(CMAKE_FIND_PACKAGE_PREFER_CONFIG TRUE)

    # https://stackoverflow.com/a/71714947/148668
    set(CGAL_DATA_DIR "unspecified")
    find_package(CGAL CONFIG COMPONENTS Core PATHS ${cgal_SOURCE_DIR} NO_DEFAULT_PATH)
endfunction()

cgal_import_target()

if(TARGET CGAL::CGAL)
    return()
endif()

message(STATUS "Third-party: creating target 'CGAL::CGAL'")

include(FetchContent)
FetchContent_Declare(
    cgal
    URL https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.12.2/CGAL-4.12.2.tar.xz
    URL_MD5 c94a0081c3836fd01ccb4d1e8bdd5d4f
    # URL https://github.com/CGAL/cgal/releases/download/v5.2.1/CGAL-5.2.1-library.tar.xz
    # URL_MD5 c1c3a9abe9106b5f3ff8dccaf2ddc0b7
)
FetchContent_GetProperties(cgal)
if(cgal_POPULATED)
    return()
endif()
FetchContent_Populate(cgal)

function(cgal_import_target)
    macro(ignore_package NAME)
        include(CMakePackageConfigHelpers)
        file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${NAME}/${NAME}Config.cmake "")
        write_basic_package_version_file(
            ${CMAKE_CURRENT_BINARY_DIR}/${NAME}/${NAME}ConfigVersion.cmake
            VERSION 1.71.0
            COMPATIBILITY AnyNewerVersion
        )
        set(${NAME}_DIR ${CMAKE_CURRENT_BINARY_DIR}/${NAME} CACHE PATH "")
        set(${NAME}_ROOT ${CMAKE_CURRENT_BINARY_DIR}/${NAME} CACHE PATH "")
    endmacro()

    # Prefer Config mode before Module mode to prevent CGAL from loading its own FindBoost.cmake
    set(CMAKE_FIND_PACKAGE_PREFER_CONFIG TRUE)

    include(boost)
    ignore_package(Boost)

    find_package(CGAL CONFIG COMPONENTS Core PATHS ${cgal_SOURCE_DIR} NO_DEFAULT_PATH)
endfunction()

cgal_import_target()

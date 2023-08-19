if(TARGET gmp::gmp)
    return()
endif()

# Download precompiled .dll on Windows
if(WIN32)
  include(gmp_mpfr)
  # Find_package will look for our downloaded lib on Windows, and system-wide on Linux/macOS
  find_package(GMP REQUIRED)
else()
  message(STATUS "Third-party: creating target 'gmp::gmp'")

  # SERIOUSLY !?! CMAKE and configure use transposed definitions of "build" and
  # "host"?
  #
  # https://cmake.org/cmake/help/latest/variable/CMAKE_SYSTEM_NAME.html#variable:CMAKE_SYSTEM_NAME
  # https://gcc.gnu.org/onlinedocs/gccint/Configure-Terms.html
  #
  # Seems these aren't to be trusted much
  # https://gitlab.kitware.com/cmake/cmake/-/issues/20989
  if(APPLE)
    # https://gmplib.org/list-archives/gmp-discuss/2020-November/006607.html
    if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" AND CMAKE_OSX_ARCHITECTURES STREQUAL "arm64")
      set(gmp_BUILD "x86_64-apple-darwin")
      set(gmp_HOST "arm64-apple-darwin")
      set(gmp_CFLAGS "--target=arm64-apple-darwin")
      set(gmp_LDFLAGS "-arch arm64")
      message(STATUS "GMP Recipe notices building on ${gmp_BUILD} for ${gmp_HOST}")
    elseif(CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64" AND CMAKE_OSX_ARCHITECTURES STREQUAL "x86_64")
      set(gmp_HOST "x86_64-apple-darwin")
      set(gmp_BUILD "arm64-apple-darwin")
      set(gmp_CFLAGS "--target=x86_64-apple-darwin13.0.0")
      set(gmp_LDFLAGS "")
      message(STATUS "GMP Recipe notices building on ${gmp_HOST} for ${gmp_BUILD}")
    endif()
  endif()

  include(FetchContent)
  include(ProcessorCount)
  ProcessorCount(Ncpu)
  include(ExternalProject)
  set(prefix ${FETCHCONTENT_BASE_DIR}/gmp)
  set(gmp_INSTALL ${prefix}/install)
  set(gmp_LIB_DIR ${gmp_INSTALL}/lib)
  set(gmp_LIBRARY 
    ${gmp_LIB_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}gmp${CMAKE_STATIC_LIBRARY_SUFFIX}
    ${gmp_LIB_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}gmpxx${CMAKE_STATIC_LIBRARY_SUFFIX}
    )
  set(gmp_INCLUDE_DIR ${gmp_INSTALL}/include)

  # Try to use CONFIGURE_HANDLED_BY_BUILD ON to avoid constantly reconfiguring
  if(${CMAKE_VERSION} VERSION_LESS 3.20)
    # CMake < 3.20, do not use any extra option
    set(gmp_ExternalProject_Add_extra_options)
  else()
    # CMake >= 3.20
    set(gmp_ExternalProject_Add_extra_options "CONFIGURE_HANDLED_BY_BUILD;ON")
  endif()

  ExternalProject_Add(gmp
    PREFIX ${prefix}
    URL  https://github.com/alisw/GMP/archive/refs/tags/v6.2.1.tar.gz
    URL_MD5 f060ad4e762ae550d16f1bb477aadba5
    UPDATE_DISCONNECTED true  # need this to avoid constant rebuild
    PATCH_COMMAND 
      curl "https://gist.githubusercontent.com/alecjacobson/d34d9307c17d1b853571699b9786e9d1/raw/8d14fc21cb7654f51c2e8df4deb0f82f9d0e8355/gmp-patch" "|" git apply -v
    ${gmp_ExternalProject_Add_extra_options}
    CONFIGURE_COMMAND 
     ${CMAKE_COMMAND} -E env
      CFLAGS=${gmp_CFLAGS}
     LDFLAGS=${gmp_LDFLAGS}
      ${prefix}/src/gmp/configure 
      --disable-debug --disable-dependency-tracking --enable-cxx --with-pic
      --prefix=${gmp_INSTALL}
      --build=${gmp_BUILD}
      --host=${gmp_HOST}
      --disable-shared
    BUILD_COMMAND make -j${Ncpu}
    INSTALL_COMMAND make -j${Ncpu} install
    INSTALL_DIR ${gmp_INSTALL}
    TEST_COMMAND ""
    BUILD_BYPRODUCTS ${gmp_LIBRARY}
  )
  ExternalProject_Get_Property(gmp SOURCE_DIR)
  set(gmp_LIBRARIES ${gmp_LIBRARY})
  add_library(gmp::gmp INTERFACE IMPORTED GLOBAL)
  file(MAKE_DIRECTORY ${gmp_INCLUDE_DIR})  # avoid race condition
  target_include_directories(gmp::gmp INTERFACE ${gmp_INCLUDE_DIR})
  target_link_libraries(gmp::gmp INTERFACE "${gmp_LIBRARIES}")  # need the quotes to expand list
  add_dependencies(gmp::gmp gmp)

endif()

if(NOT TARGET gmp::gmp)
    message(FATAL_ERROR "Creation of target 'gmp::gmp' failed")
endif()

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
    URL  https://gmplib.org/download/gmp/gmp-6.2.1.tar.xz
    URL_MD5 0b82665c4a92fd2ade7440c13fcaa42b
    UPDATE_DISCONNECTED true  # need this to avoid constant rebuild
    PATCH_COMMAND 
      curl "https://gmplib.org/repo/gmp/raw-rev/5f32dbc41afc" "|" git apply -v
    ${gmp_ExternalProject_Add_extra_options}
    CONFIGURE_COMMAND 
      ${prefix}/src/gmp/configure 
      --disable-debug --disable-dependency-tracking --enable-cxx --with-pic
      --prefix=${gmp_INSTALL}
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

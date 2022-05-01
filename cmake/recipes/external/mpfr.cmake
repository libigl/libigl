# Expects
#   gmp_INCLUDE_DIR
#   gmp_LIB_DIR
#   gmp_LIBRARIES
if(TARGET mpfr::mpfr)
    return()
endif()

# Download precompiled .dll on Windows
if(WIN32)
  include(gmp_mpfr)
  # Find_package will look for our downloaded lib on Windows, and system-wide on Linux/macOS
  find_package(MPFR REQUIRED)
else()
  message(STATUS "Third-party: creating target 'mpfr::mpfr'")

  include(FetchContent)
  include(ProcessorCount)
  ProcessorCount(Ncpu)
  include(ExternalProject)
  set(prefix ${FETCHCONTENT_BASE_DIR}/mpfr)
  set(mpfr_INSTALL ${prefix}/install)
  set(mpfr_LIBRARY ${mpfr_INSTALL}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}mpfr${CMAKE_STATIC_LIBRARY_SUFFIX})
  set(mpfr_INCLUDE_DIR ${mpfr_INSTALL}/include)

  # Try to use CONFIGURE_HANDLED_BY_BUILD ON to avoid constantly reconfiguring
  if(${CMAKE_VERSION} VERSION_LESS 3.20)
    # CMake < 3.20, do not use any extra option
    set(mpfr_ExternalProject_Add_extra_options)
  else()
    # CMake >= 3.20
    set(mpfr_ExternalProject_Add_extra_options "CONFIGURE_HANDLED_BY_BUILD;ON")
  endif()

  ExternalProject_Add(mpfr
    PREFIX ${prefix}
    DEPENDS gmp
    URL  https://ftp.gnu.org/gnu/mpfr/mpfr-4.1.0.tar.xz
    URL_MD5 bdd3d5efba9c17da8d83a35ec552baef
    UPDATE_DISCONNECTED true  # need this to avoid constant rebuild
    ${mpfr_ExternalProject_Add_extra_options} # avoid constant reconfigure
    CONFIGURE_COMMAND 
      ${prefix}/src/mpfr/configure 
      --disable-debug --disable-dependency-tracking  --disable-silent-rules --enable-cxx --with-pic
      --with-gmp-include=${gmp_INCLUDE_DIR} --with-gmp-lib=${gmp_LIB_DIR}
      --disable-shared
      --prefix=${mpfr_INSTALL}
      --disable-shared
    BUILD_COMMAND make -j${Ncpu}
    INSTALL_COMMAND make -j${Ncpu} install
    INSTALL_DIR ${mpfr_INSTALL}
    TEST_COMMAND ""
    BUILD_BYPRODUCTS ${mpfr_LIBRARY}
  )
  #PATCH_COMMAND  curl "https://raw.githubusercontent.com/Homebrew/formula-patches/03cf8088210822aa2c1ab544ed58ea04c897d9c4/libtool/configure-big_sur.diff" "|" sed -e "s/configure.orig/configure/g" "|" git apply -v
  ExternalProject_Get_Property(mpfr SOURCE_DIR)
  set(mpfr_LIBRARIES ${mpfr_LIBRARY})
  add_library(mpfr::mpfr INTERFACE IMPORTED GLOBAL)
  file(MAKE_DIRECTORY ${mpfr_INCLUDE_DIR})  # avoid race condition
  target_include_directories(mpfr::mpfr INTERFACE ${mpfr_INCLUDE_DIR})
  target_link_libraries(mpfr::mpfr INTERFACE "${mpfr_LIBRARIES}")  # need the quotes to expand list
  # This is necessary to ensure that mpfr appears before gmp in link order.
  # Otherwise undefined reference errors occur at link time on Linux with gcc
  target_link_libraries(mpfr::mpfr INTERFACE "${gmp_LIBRARIES}") 
  add_dependencies(mpfr::mpfr mpfr)
endif()

if(NOT TARGET mpfr::mpfr)
    message(FATAL_ERROR "Creation of target 'mpfr::mpfr' failed")
endif()

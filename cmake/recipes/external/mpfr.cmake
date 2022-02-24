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

  set(mpfr_INSTALL ${PROJECT_BINARY_DIR}/mpfr-install)
  set(mpfr_LIBRARY ${mpfr_INSTALL}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}mpfr${CMAKE_STATIC_LIBRARY_SUFFIX})
  set(mpfr_INCLUDE_DIR ${mpfr_INSTALL}/include)
  include(ProcessorCount)
  ProcessorCount(Ncpu)
  include(ExternalProject)
  ExternalProject_Add(mpfr
  DEPENDS gmp
  URL  https://ftp.gnu.org/gnu/mpfr/mpfr-4.1.0.tar.xz
  URL_MD5 bdd3d5efba9c17da8d83a35ec552baef
  UPDATE_DISCONNECTED true  # need this to avoid constant rebuild
  CONFIGURE_HANDLED_BY_BUILD ON  # avoid constant reconfigure
  CONFIGURE_COMMAND 
    ${PROJECT_BINARY_DIR}/mpfr-prefix/src/mpfr/configure 
    --disable-debug --disable-dependency-tracking  --disable-silent-rules --enable-cxx --with-pic
    --with-gmp-include=${gmp_INCLUDE_DIR} --with-gmp-lib=${gmp_LIB_DIR}
    --disable-shared
    --prefix=${mpfr_INSTALL}
    --disable-shared
  BUILD_COMMAND make -j${Ncpu}
  INSTALL_COMMAND make -j${Ncpu} install
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
  add_dependencies(mpfr::mpfr mpfr)
endif()

if(NOT TARGET mpfr::mpfr)
    message(FATAL_ERROR "Creation of target 'mpfr::mpfr' failed")
endif()

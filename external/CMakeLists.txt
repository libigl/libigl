cmake_minimum_required(VERSION 3.1)
project(libigl_external)

include(ExternalProject)

if(WIN32)
  # Note: It IS important to download different files on different OS's:
  # on Unix-like systems, we need the file persmissions (only available in the .tar.gz),
  # while on Windows, we need CR/LF line feeds (only available in the .zip)

  # Boost 1.65 is problematic because of https://github.com/boostorg/build/issues/230
  # Boost 1.66 is problematic because of https://github.com/boostorg/build/issues/273
  if(UNIX)
    set(BOOST_URL "https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.gz")
    set(BOOST_HASH SHA256=a13de2c8fbad635e6ba9c8f8714a0e6b4264b60a29b964b940a22554705b6b60)
    set(BOOST_BOOTSTRAP_COMMAND ./bootstrap.sh)
    set(BOOST_B2_COMMAND ./b2)
  elseif(WIN32)
    set(BOOST_URL "https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.zip")
    set(BOOST_HASH SHA256=d1775aef807bd6a14077b450cce2950d8eacd86aaf51015f79e712917f8cc3c2)
    set(BOOST_BOOTSTRAP_COMMAND cmd /C bootstrap.bat msvc)
    set(BOOST_B2_COMMAND b2.exe)
  endif()

  # We only support x64 builds
  set(BOOST_ADDRESS_MODEL 64)

  if(MSVC_VERSION EQUAL 1500) #VS2008
    set(BOOST_TOOLSET "--toolset=msvc-9.0")
  elseif(MSVC_VERSION EQUAL 1600) #VS2010
    set(BOOST_TOOLSET "--toolset=msvc-10.0")
  elseif(MSVC_VERSION EQUAL 1700) #VS2012
    set(BOOST_TOOLSET "--toolset=msvc-11.0")
  elseif(MSVC_VERSION EQUAL 1800) #VS2013
    set(BOOST_TOOLSET "--toolset=msvc-12.0")
  elseif(MSVC_VERSION EQUAL 1900) #VS2015
    set(BOOST_TOOLSET "--toolset=msvc-14.0")
  elseif(NOT (MSVC_VERSION LESS 1910)) #VS2017
    set(BOOST_TOOLSET "--toolset=msvc-14.1")
  endif()

  # Boost
  set(BOOST_DIR "${CMAKE_CURRENT_LIST_DIR}/boost")
  ExternalProject_Add(Boost_Project
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND ""

    PREFIX "${BOOST_DIR}"
    BUILD_IN_SOURCE 1

    URL ${BOOST_URL}
    URL_HASH ${BOOST_HASH}

    UPDATE_COMMAND ""
    PATCH_COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_LIST_DIR}/os.jam" "tools/build/src/util/os.jam"
    CONFIGURE_COMMAND ${BOOST_BOOTSTRAP_COMMAND} --prefix=${BOOST_DIR}/lib
    BUILD_COMMAND ${BOOST_B2_COMMAND} install -j8 -d0 --prefix=${BOOST_DIR} --with-filesystem --with-system --with-thread --with-serialization --with-date_time address-model=${BOOST_ADDRESS_MODEL} link=static runtime-link=static ${BOOST_TOOLSET}
    INSTALL_COMMAND ""
  )
endif()

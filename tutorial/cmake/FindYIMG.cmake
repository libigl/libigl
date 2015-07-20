# - Try to find the YIMG library
# Once done this will define
#
#  YIMG_FOUND - system has YIMG
#  YIMG_INCLUDE_DIR - the YIMG include directory
#  YIMG_SOURCES - the YIMG source files

FIND_PATH(YIMG_INCLUDE_DIR YImage.hpp
   /usr/include
   /usr/local/include
   /opt/local/include
   $ENV{LIBIGL}/external/yimg
   ../external/yimg/
   ../../external/yimg/
   ../../external/yimg/
   ../libigl/external/yimg/
   ../../libigl/external/yimg/
   ../../../libigl/external/yimg/

)

set(YIMG_SOURCES
${YIMG_INCLUDE_DIR}/YImage.cpp)

if(YIMG_INCLUDE_DIR)
   message(STATUS "Found YIMG: ${YIMG_INCLUDE_DIR}")
   set(YIMG_FOUND "YES")
else()
  set(YIMG_FOUND "NO")
  if (NOT YIMG_FIND_QUIETLY)
      message(FATAL_ERROR "could NOT find YIMG")
  endif (NOT YIMG_FIND_QUIETLY)
endif()

MARK_AS_ADVANCED(YIMG_INCLUDE_DIR YIMG_LIBRARIES)

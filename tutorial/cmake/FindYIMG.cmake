# - Try to find the YIMG library
# Once done this will define
#
#  YIMG_FOUND - system has YIMG
#  YIMG_INCLUDE_DIR - the YIMG include directory
#  YIMG_LIBRARIES - the YIMG libraries

FIND_PATH(YIMG_INCLUDE_DIR YImage.hpp
   /usr/include
   /usr/local/include
   /opt/local/include
   $ENV{LIBIGL}/external/yimg
)

set(YIMG_LIB_DIRS
   /usr/include
   /usr/local/include
   /opt/local/include
   $ENV{LIBIGL}/external/yimg)
FIND_LIBRARY( YIMG_LIBRARIES NAMES yimg PATHS ${YIMG_LIB_DIRS})

if(YIMG_INCLUDE_DIR AND YIMG_LIBRARIES)
   message(STATUS "Found YIMG: ${YIMG_INCLUDE_DIR}")
   set(YIMG_FOUND "YES")
else()
   message(WARNING "could NOT find YIMG")
   set(YIMG_FOUND "NO")
endif()

MARK_AS_ADVANCED(YIMG_INCLUDE_DIR YIMG_LIBRARIES)



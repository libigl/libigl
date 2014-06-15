# - Try to find the LIBIGL library
# Once done this will define
#
#  LIBIGL_FOUND - system has LIBIGL
#  LIBIGL_INCLUDE_DIR - the LIBIGL include directory
#  LIBIGL_SOURCES - the LIBIGL source files


if(LIBIGL_INCLUDE_DIR AND LIBIGL_SOURCES)
   set(LIBIGL_FOUND TRUE)
else(LIBIGL_INCLUDE_DIR AND LIBIGL_SOURCES)

FIND_PATH(LIBIGL_INCLUDE_DIR igl/readOBJ.h
   /usr/include
   /usr/local/include
   $ENV{LIBIGLROOT}/include
   $ENV{LIBIGL_ROOT}/include
   $ENV{LIBIGL_DIR}/include
   $ENV{LIBIGL_DIR}/inc
   ${PROJECT_SOURCE_DIR}/../libigl/include
   ${PROJECT_SOURCE_DIR}/../../libigl/include
   ${PROJECT_SOURCE_DIR}/../../include
)

if(LIBIGL_INCLUDE_DIR)
   set(LIBIGL_FOUND TRUE)
#   add_definitions(-DIGL_HEADER_ONLY)
   set(LIBIGL_SOURCES
      ${LIBIGL_INCLUDE_DIR}/igl/viewer/Viewer.cpp
   )
endif(LIBIGL_INCLUDE_DIR)

if(LIBIGL_FOUND)
   if(NOT LIBIGL_FIND_QUIETLY)
      message(STATUS "Found LIBIGL: ${LIBIGL_INCLUDE_DIR}")
   endif(NOT LIBIGL_FIND_QUIETLY)
else(LIBIGL_FOUND)
   if(LIBIGL_FIND_REQUIRED)
      message(FATAL_ERROR "could NOT find LIBIGL")
   endif(LIBIGL_FIND_REQUIRED)
endif(LIBIGL_FOUND)

MARK_AS_ADVANCED(LIBIGL_INCLUDE_DIR LIBIGL_LIBRARIES IGL_VIEWER_SOURCES)

endif(LIBIGL_INCLUDE_DIR AND LIBIGL_SOURCES)

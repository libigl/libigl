# - Try to find the LIBIGL library
# Once done this will define
#
#  LIBIGL_FOUND - system has LIBIGL
#  LIBIGL_INCLUDE_DIR - **the** LIBIGL include directory
#  LIBIGL_INCLUDE_DIRS - LIBIGL include directories
#  LIBIGL_SOURCES - the LIBIGL source files
if(NOT LIBIGL_FOUND)

FIND_PATH(LIBIGL_INCLUDE_DIR igl/readOBJ.h
   ${PROJECT_SOURCE_DIR}/../../include
   ${PROJECT_SOURCE_DIR}/../include
   ${PROJECT_SOURCE_DIR}/include
   ${PROJECT_SOURCE_DIR}/../external/libigl/include
   ${PROJECT_SOURCE_DIR}/../../external/libigl/include
   $ENV{LIBIGL}/include
   $ENV{LIBIGLROOT}/include
   $ENV{LIBIGL_ROOT}/include
   $ENV{LIBIGL_DIR}/include
   $ENV{LIBIGL_DIR}/inc
   /usr/include
   /usr/local/include
   /usr/local/igl/libigl/include
)


if(LIBIGL_INCLUDE_DIR)
   set(LIBIGL_FOUND TRUE)
   set(LIBIGL_INCLUDE_DIRS ${LIBIGL_INCLUDE_DIR}  ${LIBIGL_INCLUDE_DIR}/../external/Singular_Value_Decomposition)
   #set(LIBIGL_SOURCES
   #   ${LIBIGL_INCLUDE_DIR}/igl/viewer/Viewer.cpp
   #)
endif()

endif()

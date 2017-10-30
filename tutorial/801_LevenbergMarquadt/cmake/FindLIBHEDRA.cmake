# - Try to find the LIBHEDRA library
# Once done this will define
#
#  LIBHEDRA_FOUND - system has LIBHEDRA
#  LIBHEDRA_INCLUDE_DIR - **the** LIBHEDRA include directory
#  LIBHEDRA_INCLUDE_DIRS - LIBHEDRA include directories
#  LIBHEDRAL_SOURCES - the LIBHEDRA source files
if(NOT LIBHEDRA_FOUND)
message("hello")

FIND_PATH(LIBHEDRA_INCLUDE_DIR hedra/polygonal_read_OFF.h
   ${PROJECT_SOURCE_DIR}/../../include
   ${PROJECT_SOURCE_DIR}/../include
   ${PROJECT_SOURCE_DIR}/include
   /usr/include
   /usr/local/include
)

if(LIBHEDRA_INCLUDE_DIR)
   set(LIBHEDRA_FOUND TRUE)
   set(LIBHEDRA_INCLUDE_DIRS ${LIBHEDRA_INCLUDE_DIR})
endif()

endif()

# - Try to find the TINYXML2 library
# Once done this will define
#
#  TINYXML2_FOUND - system has TINYXML2
#  TINYXML2_INCLUDE_DIR - the TINYXML2 include directory
#  TINYXML2_SOURCES - the TINYXML2 source files

FIND_PATH(TINYXML2_INCLUDE_DIR tinyxml2.h
   /usr/include
   /usr/local/include
   ${PROJECT_SOURCE_DIR}/../libigl/external/tinyxml2/
   ${PROJECT_SOURCE_DIR}/../../external/tinyxml2/
)

set(TINYXML2_SOURCES ${TINYXML2_INCLUDE_DIR}/tinyxml2.cpp)

if(TINYXML2_INCLUDE_DIR)
   message(STATUS "Found TINYXML2: ${TINYXML2_INCLUDE_DIR}")
else(TINYXML2_INCLUDE_DIR)
   message(FATAL_ERROR "could NOT find TINYXML2")
endif(TINYXML2_INCLUDE_DIR)

MARK_AS_ADVANCED(TINYXML2_INCLUDE_DIR TINYXML2_LIBRARIES TINYXML2_SOURCES)

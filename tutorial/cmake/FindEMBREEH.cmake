#
# Try to find EMBREE header files
# Once done this will define
#
# EMBREE_FOUND           - system has EMBREE
# EMBREE_INCLUDE_DIRS    - the EMBREE include directories

FIND_PATH(EMBREE_INCLUDE_DIR embree2/rtcore.h
	  PATHS
		${PROJECT_SOURCE_DIR}/../../external/embree/include
		${PROJECT_SOURCE_DIR}/../external/embree/include
		${PROJECT_SOURCE_DIR}/../libigl/external/embree/include
		NO_DEFAULT_PATH
    )

if(EMBREE_INCLUDE_DIR)
	set(EMBREE_FOUND TRUE)
endif(EMBREE_INCLUDE_DIR)

IF (EMBREE_FOUND)
   message(STATUS "Found EMBREE: ${EMBREE_INCLUDE_DIR}")

   SET(EMBREE_INCLUDE_DIRS ${EMBREE_INCLUDE_DIR} ${EMBREE_INCLUDE_DIR}/embree)
ELSE (EMBREE_FOUND)
    message(STATUS "could NOT find EMBREE")
ENDIF (EMBREE_FOUND)

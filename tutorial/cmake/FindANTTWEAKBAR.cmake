#
# Try to find AntTweakBar library and include path.
# Once done this will define
#
# ANT_TWEAK_BAR_FOUND
# ANT_TWEAK_BAR_INCLUDE_DIR
# ANT_TWEAK_BAR_LIBRARY
#

IF (WIN32)
	IF( CMAKE_SIZEOF_VOID_P EQUAL 8 )
	  SET( BITS "64" )
    ELSE( CMAKE_SIZEOF_VOID_P EQUAL 8 )
      SET( BITS "" )
    ENDIF( CMAKE_SIZEOF_VOID_P EQUAL 8 )

	FIND_PATH( ANT_TWEAK_BAR_INCLUDE_DIR AntTweakBar.h
      PATHS
		${PROJECT_SOURCE_DIR}/../../external/AntTweakBar/include
		${PROJECT_SOURCE_DIR}/../external/AntTweakBar/include
		$ENV{ANT_TWEAK_BAR_ROOT}/include
		DOC "The directory where AntTweakBar.h resides")

    FIND_LIBRARY( ANT_TWEAK_BAR_LIBRARY AntTweakBar${BITS}
        PATHS
		${PROJECT_SOURCE_DIR}/../../external/AntTweakBar/lib
		${PROJECT_SOURCE_DIR}/../external/AntTweakBar/lib
                $ENV{ANT_TWEAK_BAR_ROOT}/lib
                DOC "The AntTweakBar library")
ELSE (WIN32)

FIND_PATH(ANT_TWEAK_BAR_INCLUDE_DIR AntTweakBar.h
      PATHS
	    ${LIBIGL_INCLUDE_DIR}/../external/AntTweakBar/include/
      ${PROJECT_SOURCE_DIR}/../../external/AntTweakBar/include/
      ${PROJECT_SOURCE_DIR}/../external/AntTweakBar/include/
      /usr/local/include
      /usr/X11/include
      /usr/include
      NO_DEFAULT_PATH)

FIND_LIBRARY( ANT_TWEAK_BAR_LIBRARY AntTweakBar
  PATHS
		${LIBIGL_INCLUDE_DIR}/../external/AntTweakBar/lib
    ${PROJECT_SOURCE_DIR}/../../external/AntTweakBar/lib
    ${PROJECT_SOURCE_DIR}/../external/AntTweakBar/lib
    /usr/local
    /usr/X11
    /usr
  PATH_SUFFIXES
    a
    lib64
    lib
    dylib
    NO_DEFAULT_PATH
)

ENDIF (WIN32)

set(ANT_TWEAK_BAR_INCLUDE_DIR ${ANT_TWEAK_BAR_INCLUDE_DIR} ${ANT_TWEAK_BAR_INCLUDE_DIR}/../src/)

# message(FATAL_ERROR ${ANT_TWEAK_BAR_LIBRARY})

if(ANT_TWEAK_BAR_INCLUDE_DIR AND ANT_TWEAK_BAR_LIBRARY)
	message(STATUS "Found ANTTWEAKBAR: ${ANT_TWEAK_BAR_INCLUDE_DIR}")
else(ANT_TWEAK_BAR_INCLUDE_DIR AND ANT_TWEAK_BAR_LIBRARY)
	message(FATAL_ERROR "could NOT find ANTTWEAKBAR")
endif(ANT_TWEAK_BAR_INCLUDE_DIR AND ANT_TWEAK_BAR_LIBRARY)

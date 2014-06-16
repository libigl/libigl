#
# Try to find AntTweakBar library and include path.
# Once done this will define
#
# ANT_TWEAK_BAR_FOUND
# ANT_TWEAK_BAR_INCLUDE_DIR
# ANT_TWEAK_BAR_LIBRARY
#

FIND_PATH(ANT_TWEAK_BAR_INCLUDE_DIR AntTweakBar.h
      PATHS
	    ${PROJECT_SOURCE_DIR}/../libigl/external/AntTweakBar/include/
      ${PROJECT_SOURCE_DIR}/../../external/AntTweakBar/include/
      /usr/local/include
      /usr/X11/include
      /usr/include
      NO_DEFAULT_PATH)

set(ANT_TWEAK_BAR_INCLUDE_DIR ${ANT_TWEAK_BAR_INCLUDE_DIR} ${ANT_TWEAK_BAR_INCLUDE_DIR}/../src/)

FIND_LIBRARY( ANT_TWEAK_BAR_LIBRARY AntTweakBar
  PATHS
		${PROJECT_SOURCE_DIR}/../libigl/external/AntTweakBar/lib
    ${PROJECT_SOURCE_DIR}/../../external/AntTweakBar/lib
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

# message(FATAL_ERROR ${ANT_TWEAK_BAR_LIBRARY})

if(ANT_TWEAK_BAR_INCLUDE_DIR AND ANT_TWEAK_BAR_LIBRARY)
	message(STATUS "Found ANTTWEAKBAR: ${ANT_TWEAK_BAR_INCLUDE_DIR}")
else(ANT_TWEAK_BAR_INCLUDE_DIR AND ANT_TWEAK_BAR_LIBRARY)
	message(FATAL_ERROR "could NOT find ANTTWEAKBAR")
endif(ANT_TWEAK_BAR_INCLUDE_DIR AND ANT_TWEAK_BAR_LIBRARY)

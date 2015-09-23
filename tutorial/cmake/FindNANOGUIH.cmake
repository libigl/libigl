#
# Try to find NANOGUI library and include path.
# Once done this will define
#
# NANOGUI_FOUND
# NANOGUI_INCLUDE_DIR
# NANOGUI_LIBRARY
#

if(NOT NANOGUI_FOUND)

FIND_PATH(NANOGUI_INCLUDE_DIR nanogui/nanogui.h
  PATHS
    ${PROJECT_SOURCE_DIR}/../../external/nanogui/include
    ${PROJECT_SOURCE_DIR}/../external/nanogui/include
    ${PROJECT_SOURCE_DIR}/external/nanogui/include
    ${PROJECT_SOURCE_DIR}/../../libigl/external/nanogui/include
    ${PROJECT_SOURCE_DIR}/../libigl/external/nanogui/include
    ${PROJECT_SOURCE_DIR}/libigl/external/nanogui/include
    /usr/local/include
    /usr/X11/include
    /usr/include
    /opt/local/include
    NO_DEFAULT_PATH
    )

SET(NANOGUI_FOUND "NO")
IF (NANOGUI_INCLUDE_DIR)
	SET(NANOGUI_FOUND "YES")
  SET(NANOGUI_INCLUDE_DIRS
         ${NANOGUI_INCLUDE_DIR}
         ${NANOGUI_INCLUDE_DIR}/../ext/nanovg/src
         ${NANOGUI_INCLUDE_DIR}/../ext/glfw/include
         )

ENDIF (NANOGUI_INCLUDE_DIR)

if(NANOGUI_FOUND)
  message(STATUS "Found NANOGUI: ${NANOGUI_INCLUDE_DIR} -- HEADERS ONLY")
else(NANOGUI_FOUND)
  if (NOT NANOGUIH_FIND_QUIETLY)
    message(FATAL_ERROR "could NOT find NANOGUI")
  endif (NOT NANOGUIH_FIND_QUIETLY)
endif(NANOGUI_FOUND)

endif(NOT NANOGUI_FOUND)

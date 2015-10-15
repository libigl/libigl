#
# Try to find GLFW library and include path.
# Once done this will define
#
# GLFW_FOUND
# GLFW_INCLUDE_DIR
# GLFW_LIBRARIES
#

if(NOT GLFW_FOUND)

FIND_PATH(GLFW_INCLUDE_DIR GLFW/glfw3.h
  PATHS
    ${PROJECT_SOURCE_DIR}/../../external/glfw/include
    ${PROJECT_SOURCE_DIR}/../external/glfw/include
    ${PROJECT_SOURCE_DIR}/external/glfw/include
    ${PROJECT_SOURCE_DIR}/../../libigl/external/glfw/include
    ${PROJECT_SOURCE_DIR}/../libigl/external/glfw/include
    ${PROJECT_SOURCE_DIR}/libigl/external/glfw/include
    ${PROJECT_SOURCE_DIR}/../../../libigl/external/nanogui/ext/glfw/include
    ${PROJECT_SOURCE_DIR}/../../libigl/external/nanogui/ext/glfw/include
    ${PROJECT_SOURCE_DIR}/../libigl/external/nanogui/ext/glfw/include
    ${PROJECT_SOURCE_DIR}/libigl/external/nanogui/ext/glfw/include
    /usr/local/include
    /usr/X11/include
    /usr/include
    /opt/local/include
    NO_DEFAULT_PATH
    )

SET(GLFW_FOUND "NO")
IF (GLFW_INCLUDE_DIR)
	SET(GLFW_FOUND "YES")
ENDIF (GLFW_INCLUDE_DIR)

if(GLFW_FOUND)
  message(STATUS "Found GLFW: ${GLFW_INCLUDE_DIR} -- HEADERS ONLY")
else(GLFW_FOUND)
  if (NOT GLFW_FIND_QUIETLY)
    message(WARNING "could NOT find GLFW")
  endif (NOT GLFW_FIND_QUIETLY)
endif(GLFW_FOUND)

endif(NOT GLFW_FOUND)

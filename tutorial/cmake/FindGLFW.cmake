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
    ${PROJECT_SOURCE_DIR}/../../libigl/external/nanogui/ext/glfw/include
    ${PROJECT_SOURCE_DIR}/../libigl/external/nanogui/ext/glfw/include
    ${PROJECT_SOURCE_DIR}/libigl/external/nanogui/ext/glfw/include
    /usr/local/include
    /usr/X11/include
    /usr/include
    /opt/local/include
    NO_DEFAULT_PATH
    )

FIND_LIBRARY( GLFW_LIBRARIES NAMES glfw glfw3
  PATHS
    ${PROJECT_SOURCE_DIR}/../../external/glfw/src
    ${PROJECT_SOURCE_DIR}/../external/glfw/src
    ${PROJECT_SOURCE_DIR}/external/glfw/src
    ${PROJECT_SOURCE_DIR}/../../libigl/external/glfw/src
    ${PROJECT_SOURCE_DIR}/../libigl/external/glfw/src
    ${PROJECT_SOURCE_DIR}/libigl/external/glfw/src
    ${PROJECT_SOURCE_DIR}/../../external/glfw/lib/x64
    ${PROJECT_SOURCE_DIR}/../external/glfw/lib/x64
    ${PROJECT_SOURCE_DIR}/external/glfw/lib/x64
    ${PROJECT_SOURCE_DIR}/../../libigl/external/glfw/lib/x64
    ${PROJECT_SOURCE_DIR}/../libigl/external/glfw/lib/x64
    ${PROJECT_SOURCE_DIR}/libigl/external/glfw/lib/x64
    /usr/local
    /usr/X11
    /usr
    PATH_SUFFIXES
    a
    lib64
    lib
    NO_DEFAULT_PATH
)

SET(GLFW_FOUND "NO")
IF (GLFW_INCLUDE_DIR AND GLFW_LIBRARIES)
	SET(GLFW_FOUND "YES")
ENDIF (GLFW_INCLUDE_DIR AND GLFW_LIBRARIES)

if(GLFW_FOUND)
  message(STATUS "Found GLFW: ${GLFW_INCLUDE_DIR}")
else(GLFW_FOUND)
  if (NOT GLFW_FIND_QUIETLY)
    message(FATAL_ERROR "could NOT find GLFW")
  endif (NOT GLFW_FIND_QUIETLY)
endif(GLFW_FOUND)

endif(NOT GLFW_FOUND)

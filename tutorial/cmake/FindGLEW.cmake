# - Try to find the GLEW library
# Once done this will define
#
#  GLEW_FOUND - system has GLEW
#  GLEW_INCLUDE_DIR - the GLEW include directory
#  GLEW_SOURCES - the GLEW source file list

FIND_PATH(GLEW_INCLUDE_DIR GL/glew.h
   ${PROJECT_SOURCE_DIR}/../../../external/nanogui/ext/glew/include
   ${PROJECT_SOURCE_DIR}/../../external/nanogui/ext/glew/include
   ${PROJECT_SOURCE_DIR}/../external/nanogui/ext/glew/include
   ${PROJECT_SOURCE_DIR}/external/nanogui/ext/glew/include
   ${PROJECT_SOURCE_DIR}/../../../libigl/external/nanogui/ext/glew/include
   ${PROJECT_SOURCE_DIR}/../../libigl/external/nanogui/ext/glew/include
   ${PROJECT_SOURCE_DIR}/../libigl/external/nanogui/ext/glew/include
   ${PROJECT_SOURCE_DIR}/libigl/external/nanogui/ext/glew/include
   /usr/include
   /usr/local/include
   $ENV{GLEWROOT}/include
   $ENV{GLEW_ROOT}/include
   $ENV{GLEW_DIR}/include
   $ENV{GLEW_DIR}/inc
   [HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\8.0\\Setup\\VC]/PlatformSDK/Include
   NO_DEFAULT_PATH
)

if(GLEW_INCLUDE_DIR)
   set(GLEW_FOUND TRUE)
endif(GLEW_INCLUDE_DIR)

if(GLEW_FOUND)
  set(GLEW_SOURCES ${GLEW_INCLUDE_DIR}/../src/glew.c)
  message(STATUS "Found GLEW: ${GLEW_INCLUDE_DIR}")
else(GLEW_FOUND)
  message(WARNING "could NOT find glew")
endif(GLEW_FOUND)

MARK_AS_ADVANCED(GLEW_INCLUDE_DIR)

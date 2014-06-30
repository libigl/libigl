# - Try to find the LIBCOMISO library
# Once done this will define
#
#  LIBCOMISO_FOUND - system has LIBCOMISO
#  LIBCOMISO_INCLUDE_DIR - the LIBCOMISO include directory
#  LIBCOMISO_LIBRARY - the LIBCOMISO binary lib

FIND_PATH(LIBCOMISO_INCLUDE_DIR CoMISo/Solver/ConstrainedSolver.hh
   /usr/include
   /usr/local/include
   $ENV{LIBCOMISOROOT}/include
   $ENV{LIBCOMISO_ROOT}/include
   $ENV{LIBCOMISO_DIR}/include
   $ENV{LIBCOMISO_DIR}/inc
   ${PROJECT_SOURCE_DIR}/../
   ${PROJECT_SOURCE_DIR}/../../
   ${PROJECT_SOURCE_DIR}/../../../
   ${PROJECT_SOURCE_DIR}/../CoMISo/
   ${PROJECT_SOURCE_DIR}/../CoMISo/include
   ${PROJECT_SOURCE_DIR}/../../CoMISo/
   ${PROJECT_SOURCE_DIR}/../../CoMISo/include
   /Users/daniele/Dropbox/igl/MIQ/src
   /Users/olkido/Documents/igl/MIQ/src
)

#message(FATAL_ERROR "${LIBCOMISO_INCLUDE_DIR}")

FIND_LIBRARY(LIBCOMISO_LIBRARY NAMES CoMISo
  PATHS
    /usr/local
    /usr/X11
    /usr
    /
    ${PROJECT_SOURCE_DIR}/../CoMISo/
    ${PROJECT_SOURCE_DIR}/../CoMISo/build/Build/lib/CoMISo/
    ${PROJECT_SOURCE_DIR}/../../CoMISo/
    ${PROJECT_SOURCE_DIR}/../../CoMISo/build/Build/lib/CoMISo/
    ${PROJECT_SOURCE_DIR}/../../../CoMISo/
    ${PROJECT_SOURCE_DIR}/../../../CoMISo/build/Build/lib/CoMISo/
    /Users/olkido/Documents/igl/MIQ/src/CoMISo/Build
)
#message(STATUS "${LIBCOMISO_LIBRARY}")

if(LIBCOMISO_INCLUDE_DIR AND LIBCOMISO_LIBRARY)

   set(LIBCOMISO_INCLUDE_DIR ${LIBCOMISO_INCLUDE_DIR} ${LIBCOMISO_INCLUDE_DIR}/CoMISo/gmm/include)

   add_definitions(-DINCLUDE_TEMPLATES)
   message(STATUS "Found LIBCOMISO: ${LIBCOMISO_INCLUDE_DIR} ${LIBCOMISO_LIBRARY}")
   set(LIBCOMISO_FOUND TRUE)
else(LIBCOMISO_INCLUDE_DIR)
   message(STATUS "could NOT find LIBCOMISO")
endif(LIBCOMISO_INCLUDE_DIR AND LIBCOMISO_LIBRARY)

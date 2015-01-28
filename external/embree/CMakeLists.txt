## ======================================================================== ##
## Copyright 2009-2013 Intel Corporation                                    ##
##                                                                          ##
## Licensed under the Apache License, Version 2.0 (the "License");          ##
## you may not use this file except in compliance with the License.         ##
## You may obtain a copy of the License at                                  ##
##                                                                          ##
##     http://www.apache.org/licenses/LICENSE-2.0                           ##
##                                                                          ##
## Unless required by applicable law or agreed to in writing, software      ##
## distributed under the License is distributed on an "AS IS" BASIS,        ##
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. ##
## See the License for the specific language governing permissions and      ##
## limitations under the License.                                           ##
## ======================================================================== ##

PROJECT(embree)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
MARK_AS_ADVANCED(CMAKE_INSTALL_PREFIX)
MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH)
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH)
MARK_AS_ADVANCED(CLEAR CMAKE_VERBOSE_MAKEFILE)

MARK_AS_ADVANCED(CMAKE_OSX_ARCHITECTURES)
MARK_AS_ADVANCED(CMAKE_OSX_DEPLOYMENT_TARGET)
MARK_AS_ADVANCED(CMAKE_OSX_SYSROOT)
MARK_AS_ADVANCED(GLUT_cocoa_LIBRARY)

SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/common/cmake ${CMAKE_MODULE_PATH})

##############################################################
# Embree configuration
##############################################################

SET(RTCORE_EXIT_ON_ERROR OFF CACHE BOOL "Exits the application if an error occurs.")
IF (RTCORE_EXIT_ON_ERROR)
 ADD_DEFINITIONS(-D__EXIT_ON_ERROR__)
ENDIF()

SET(RTCORE_RAY_MASK OFF CACHE BOOL "Enables ray mask support.")
IF (RTCORE_RAY_MASK)
 ADD_DEFINITIONS(-D__USE_RAY_MASK__)
ENDIF()

SET(RTCORE_STAT_COUNTERS OFF CACHE BOOL "Enables statistic counters.")
IF (RTCORE_STAT_COUNTERS)
 ADD_DEFINITIONS(-D__USE_STAT_COUNTERS__)
ENDIF()

SET(RTCORE_BACKFACE_CULLING OFF CACHE BOOL "Enables backface culling.")
IF (RTCORE_BACKFACE_CULLING)
 ADD_DEFINITIONS(-D__BACKFACE_CULLING__)
ENDIF()

IF (TARGET_XEON_PHI)
  UNSET(RTCORE_INTERSECTION_FILTER CACHE)
  SET(RTCORE_INTERSECTION_FILTER OFF CACHE BOOL "Enables intersection filter callback.")
ELSE ()
  SET(RTCORE_INTERSECTION_FILTER ON CACHE BOOL "Enables intersection filter callback.")
ENDIF()

IF (RTCORE_INTERSECTION_FILTER)
 ADD_DEFINITIONS(-D__INTERSECTION_FILTER__)
ENDIF()

IF (TARGET_XEON_PHI)
  UNSET(RTCORE_BUFFER_STRIDE CACHE)
  SET(RTCORE_BUFFER_STRIDE OFF CACHE BOOL "Enables support for strides and offsets of buffers.")
ELSE (TARGET_XEON_PHI)
  SET(RTCORE_BUFFER_STRIDE ON CACHE BOOL "Enables support for strides and offsets of buffers.")
ENDIF()

IF (RTCORE_BUFFER_STRIDE)
 ADD_DEFINITIONS(-D__BUFFER_STRIDE__)
ENDIF()

IF (APPLE)
  SET(RTCORE_SPINLOCKS OFF CACHE BOOL "Use spinning locks.")
ELSE()
  SET(RTCORE_SPINLOCKS ON CACHE BOOL "Use spinning locks.")
ENDIF()

IF (RTCORE_SPINLOCKS)
 ADD_DEFINITIONS(-D__SPINLOCKS__)
ENDIF()

SET(RTCORE_TASKLOGGER OFF CACHE BOOL "Allows creating scheduling diagram of tasks.")
IF (RTCORE_TASKLOGGER)
 ADD_DEFINITIONS(-D__LOG_TASKS__)
ENDIF()

SET(RTCORE_EXPORT_ALL_SYMBOLS OFF CACHE BOOL "Lets Embree shared library export all symbols.")
IF (RTCORE_EXPORT_ALL_SYMBOLS)
 ADD_DEFINITIONS(-D__EXPORT_ALL_SYMBOLS__)
ENDIF()

SET(RTCORE_FIX_RAYS OFF CACHE BOOL "Makes traversal algorithms NAN and INF safe")
IF (RTCORE_FIX_RAYS)
 ADD_DEFINITIONS(-D__FIX_RAYS__)
ENDIF()

##############################################################
# ISA configuration
##############################################################

MACRO (TEST_COMPILER)

  STRING(REGEX REPLACE " " ";" FLAGS_SSE41_LIST ${FLAGS_SSE41})
  EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER} ${FLAGS_SSE41_LIST} -c ${CMAKE_CURRENT_SOURCE_DIR}/common/cmake/dummy.cpp
                  RESULT_VARIABLE res OUTPUT_VARIABLE out ERROR_VARIABLE err)
  IF(NOT ${res} EQUAL 0)
    UNSET(TARGET_SSE41 CACHE)
    SET(TARGET_SSE41 OFF CACHE BOOL "Generate code path for SSE41")
  ENDIF()

  STRING(REGEX REPLACE " " ";" FLAGS_AVX_LIST ${FLAGS_AVX})
  EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER} ${FLAGS_AVX_LIST} -c ${CMAKE_CURRENT_SOURCE_DIR}/common/cmake/dummy.cpp
                  RESULT_VARIABLE res OUTPUT_VARIABLE out ERROR_VARIABLE err)
  IF(NOT ${res} EQUAL 0)
    UNSET(TARGET_AVX CACHE)
    SET(TARGET_AVX   OFF CACHE BOOL "Generate code path for AVX")
  ENDIF()

  STRING(REGEX REPLACE " " ";" FLAGS_AVX2_LIST ${FLAGS_AVX2})
  EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER} ${FLAGS_AVX2_LIST} -c ${CMAKE_CURRENT_SOURCE_DIR}/common/cmake/dummy.cpp
                  RESULT_VARIABLE res OUTPUT_VARIABLE out ERROR_VARIABLE err)
  IF(NOT ${res} EQUAL 0)
    UNSET(TARGET_AVX2 CACHE)
    SET(TARGET_AVX2  OFF CACHE BOOL "Generate code path for AVX2")
  ENDIF()

ENDMACRO()

SET(TARGET_SSE2  ON CACHE BOOL "Generate code path for SSE2")
SET(TARGET_SSE41 ON CACHE BOOL "Generate code path for SSE41")
SET(TARGET_AVX   ON CACHE BOOL "Generate code path for AVX")
SET(TARGET_AVX2  ON CACHE BOOL "Generate code path for AVX2")
SET(TARGET_XEON_PHI OFF CACHE BOOL "Generate code path for Xeon Phi")

IF (TARGET_SSE2 STREQUAL "ON")
 SET(ISPC_TARGETS "sse2")
ENDIF()

IF (TARGET_SSE41 STREQUAL "ON")
  SET(ISPC_TARGETS ${ISPC_TARGETS} "sse4")
ENDIF()

IF (TARGET_AVX STREQUAL "ON")
  SET(ISPC_TARGETS ${ISPC_TARGETS} "avx")
ENDIF()

IF (TARGET_AVX2 STREQUAL "ON")
  SET(ISPC_TARGETS ${ISPC_TARGETS} "avx2")
ENDIF()

SET(ISPC_TARGETS "${ISPC_TARGETS}")
STRING(REGEX REPLACE ";" "," ISPC_TARGETS "${ISPC_TARGETS}")

##############################################################
# Compiler
##############################################################

SET(COMPILER "GCC" CACHE STRING "Set to GCC, CLANG, or ICC")

IF (COMPILER STREQUAL "GCC") 
  INCLUDE (gcc)
  TEST_COMPILER()
ELSEIF (COMPILER STREQUAL "CLANG") 
  INCLUDE (clang)
  TEST_COMPILER()
ELSEIF (COMPILER STREQUAL "ICC") 
  INCLUDE (icc)
  TEST_COMPILER()
ELSE ()
  MESSAGE(FATAL_ERROR "Unknown compiler specified: " ${COMPILER})
ENDIF ()

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the build type: Debug Release RelWithDebInfo MinSizeRel." FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

SET(CMAKE_VERBOSE_MAKEFILE false CACHE BOOL "Enables verbose mode.")

INCLUDE (ispc)

##############################################################
# Search paths
##############################################################
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/common)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES_ISPC(${PROJECT_SOURCE_DIR})
LINK_DIRECTORIES(${EXECUTABLE_OUTPUT_PATH})

##############################################################
# Output paths
##############################################################

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

##############################################################
# Directories to compile
##############################################################

ADD_SUBDIRECTORY(common)
ADD_SUBDIRECTORY(kernels)
ADD_SUBDIRECTORY(tutorials)

SET (__XEON__ 1)
ADD_SUBDIRECTORY(tests)

IF (TARGET_XEON_PHI)
  SET (__XEON__ 0)
  ADD_SUBDIRECTORY(tests tests_xeonphi)
ENDIF ()

##############################################################
# Install Headers
##############################################################

INSTALL(DIRECTORY include/embree2 DESTINATION include FILES_MATCHING PATTERN "*.h")
INSTALL(DIRECTORY include/embree2 DESTINATION include FILES_MATCHING PATTERN "*.isph")


# ##################################################################
# Create Documentation. Call via 'make doc'
# ##################################################################
FIND_PACKAGE(Doxygen)
IF(DOXYGEN_FOUND)
	CONFIGURE_FILE(
		${CMAKE_CURRENT_SOURCE_DIR}/doc/dox/doxygen_config
		${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
	ADD_CUSTOM_TARGET(doc
		${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
		COMMENT "Generating Doxygen documentation" VERBATIM
		)
ENDIF(DOXYGEN_FOUND)

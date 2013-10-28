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

IF(COMMAND cmake_policy)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)

MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
MARK_AS_ADVANCED(CMAKE_INSTALL_PREFIX)
MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH)
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH)
MARK_AS_ADVANCED(USE_STAT_COUNTERS)
MARK_AS_ADVANCED(CLEAR CMAKE_VERBOSE_MAKEFILE)

##############################################################
# configuration
##############################################################

MACRO( CONFIGURE_CPU )
  IF (${TARGET_CPU} STREQUAL "xeon")
    SET(__XEON__ ON)
    SET(__XEON_PHI__ OFF)
    SET(EXT "")
  ELSEIF (${TARGET_CPU} STREQUAL "xeon_phi")
    INCLUDE (build_icc_mic)
    SET(__XEON__ OFF)
    SET(__XEON_PHI__ ON)
    SET(EXT "_knc")
  ELSE ()
    MESSAGE(FATAL_ERROR "Unknown CPU type: " ${TARGET_CPU})
  ENDIF ()
ENDMACRO( CONFIGURE_CPU )

MACRO( CONFIGURE_SPMD )
  IF (${TARGET_SPMD} STREQUAL "sse")
    SET(__XEON__ ON)
    SET(__XEON_PHI__ OFF)
    SET (EXT "_sse")
  ELSEIF (${TARGET_SPMD} STREQUAL "avx")
    SET(__XEON__ ON)
    SET(__XEON_PHI__ OFF)
    SET (EXT "_avx")
  ELSEIF (${TARGET_SPMD} STREQUAL "knc")
    SET(__XEON__ OFF)
    SET(__XEON_PHI__ ON)
    SET (EXT "_knc")
  ELSE ()
    MESSAGE(FATAL_ERROR "Unknown SPMD target: " ${TARGET_SPMD})
  ENDIF ()
ENDMACRO( CONFIGURE_SPMD )

##############################################################
# search paths
##############################################################

LINK_DIRECTORIES(${EXECUTABLE_OUTPUT_PATH})
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR} ${CMAKE_BINARY_DIR} ${PROJECT_SOURCE_DIR}/embree) 
SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

##############################################################
# Compiler
##############################################################

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

SET(CMAKE_VERBOSE_MAKEFILE false CACHE BOOL "Enables verbose mode.")

SET(COMPILER "GCC" CACHE STRING "Set to GCC, CLANG, ICC, or ICC-MIC")

IF (COMPILER STREQUAL "GCC") 
  INCLUDE (build_gcc)
ELSEIF (COMPILER STREQUAL "CLANG") 
  INCLUDE (build_clang)
ELSEIF (COMPILER STREQUAL "ICC") 
  INCLUDE (build_icc)
ELSE ()
  MESSAGE(FATAL_ERROR "Unknown compiler specified: " ${COMPILER})
ENDIF ()

##############################################################
# output paths
##############################################################

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

##############################################################
# directories to compile
##############################################################

SET(TARGET_CPU "xeon")
ADD_SUBDIRECTORY(embree)
ADD_SUBDIRECTORY(examples)

IF (BUILD_SINGLERAY_DEVICE_KNC OR BUILD_ISPC_DEVICE_KNC)
  SET(TARGET_CPU "xeon_phi")
  ADD_SUBDIRECTORY(embree embree_knc)
ENDIF ()

IF (BUILD_ISPC_DEVICE_SSE OR BUILD_TUTORIALS_SSE)
  SET(TARGET_SPMD "sse")
  ADD_SUBDIRECTORY(embree_ispc embree_ispc_sse)
ENDIF ()

IF (BUILD_ISPC_DEVICE_AVX OR BUILD_TUTORIALS_AVX)
  SET(TARGET_SPMD "avx")
  ADD_SUBDIRECTORY(embree_ispc embree_ispc_avx)
ENDIF ()

IF (BUILD_ISPC_DEVICE_KNC)
  SET(TARGET_SPMD "knc")
  ADD_SUBDIRECTORY(embree_ispc embree_ispc_knc)
ENDIF ()

FILE(WRITE "${CMAKE_BINARY_DIR}/CMakeDefines.h" "#define CMAKE_BUILD_DIR \"${CMAKE_BINARY_DIR}\"\n")

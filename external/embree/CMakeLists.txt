## ======================================================================== ##
## Copyright 2009-2012 Intel Corporation                                    ##
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

##############################################################
# Configuration
##############################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

IF(COMMAND cmake_policy)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release." FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

SET(USE_INTEL_COMPILER 0 CACHE BOOL "Set to 1 to use the Intel Compiler")
SET(SSE_VERSION "SSE4.2" CACHE INT "SSE version to use (SSSE3,SSE4.1,SSE4.2,AVX)")
SET(CMAKE_VERBOSE_MAKEFILE false)

##############################################################
# Compiler
##############################################################

IF (USE_INTEL_COMPILER)

  IF(SSE_VERSION STREQUAL "SSSE3")
    SET(SSE_FLAGS "-xssse3")
  ELSEIF(SSE_VERSION STREQUAL "SSE4.1")
    SET(SSE_FLAGS "-xsse4.1")
  ELSEIF(SSE_VERSION STREQUAL "SSE4.2")
    SET(SSE_FLAGS "-xsse4.2")
  ELSEIF(SSE_VERSION STREQUAL "AVX")
    SET(SSE_FLAGS "-xAVX")
  ELSE()
    MESSAGE(FATAL_ERROR "SSE Version ${SSE_VERSION} not supported use SSSE3, SSE4.1, SSE4.2 or AVX" )
  ENDIF()
  
  SET(CMAKE_CXX_COMPILER "icpc")
  SET(CMAKE_C_COMPILER "icc")
  SET(CMAKE_CXX_FLAGS "-Wall -fp-model fast -fvisibility-inlines-hidden -fvisibility=hidden ${SSE_FLAGS}")
  SET(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -g -O0")
  SET(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -g -O2")
  SET(CMAKE_EXE_LINKER_FLAGS "") 

ELSE(USE_INTEL_COMPILER)

  IF(SSE_VERSION STREQUAL "SSSE3")
    SET(SSE_FLAGS "-mssse3")
  ELSEIF(SSE_VERSION STREQUAL "SSE4.1")
    SET(SSE_FLAGS "-msse4.1")
  ELSEIF(SSE_VERSION STREQUAL "SSE4.2")
    SET(SSE_FLAGS "-msse4.2")
  ELSEIF(SSE_VERSION STREQUAL "AVX")
    SET(SSE_FLAGS "-mavx")
  ELSE()
    MESSAGE(FATAL_ERROR "SSE Version ${SSE_VERSION} not supported use SSSE3, SSE4.1, SSE4.2 or AVX" )
  ENDIF()

  SET(CMAKE_CXX_COMPILER "g++")
  SET(CMAKE_C_COMPILER "gcc")
  SET(CMAKE_CXX_FLAGS "-Wall -ffast-math -fvisibility-inlines-hidden -fvisibility=hidden ${SSE_FLAGS}")
  SET(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -g -O0 -ftree-ter")
  SET(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -g -O2")
  SET(CMAKE_EXE_LINKER_FLAGS "")

ENDIF(USE_INTEL_COMPILER)

##############################################################
# search paths
##############################################################

LINK_DIRECTORIES(${EXECUTABLE_OUTPUT_PATH})
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/common)

##############################################################
# output paths
##############################################################

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

##############################################################
# directories to compile
##############################################################

ADD_SUBDIRECTORY(common)
ADD_SUBDIRECTORY(rtcore)
ADD_SUBDIRECTORY(renderer)
ADD_SUBDIRECTORY(viewer)

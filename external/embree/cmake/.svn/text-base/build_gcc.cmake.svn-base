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

SET(COMPILER_TARGET "SSE4.2" CACHE INT "SIMD version to use (SSSE3,SSE4.1,SSE4.2,AVX,AVX-I,AVX2)")

IF(COMPILER_TARGET STREQUAL "SSSE3")
  SET(__SSE__ 1)
  SET(__AVX__ 0)
  SET(SIMD_FLAGS "-mssse3")
ELSEIF(COMPILER_TARGET STREQUAL "SSE4.1")
  SET(__SSE__ 1)
  SET(__AVX__ 0)
  SET(SIMD_FLAGS "-msse4.1")
ELSEIF(COMPILER_TARGET STREQUAL "SSE4.2")
  SET(__SSE__ 1)
  SET(__AVX__ 0)
  SET(SIMD_FLAGS "-msse4.2")
ELSEIF(COMPILER_TARGET STREQUAL "AVX")
  SET(__SSE__ 1)
  SET(__AVX__ 1)
  SET(SIMD_FLAGS "-mavx")
ELSEIF(COMPILER_TARGET STREQUAL "AVX-I")
  SET(__SSE__ 1)
  SET(__AVX__ 1)
  SET(SIMD_FLAGS "-mavx2")
ELSEIF(COMPILER_TARGET STREQUAL "AVX2")
  SET(__SSE__ 1)
  SET(__AVX__ 1)
  SET(SIMD_FLAGS "-mavx2")
ELSE()
  MESSAGE(FATAL_ERROR "SIMD Version ${COMPILER_TARGET} not supported use SSSE3, SSE4.1, SSE4.2, AVX, AVX-I, or AVX2" )
ENDIF()

SET(CMAKE_CXX_COMPILER "g++")
SET(CMAKE_C_COMPILER "gcc")
SET(CMAKE_CXX_FLAGS "-fPIC -ffast-math -fvisibility-inlines-hidden -fvisibility=hidden ${SIMD_FLAGS}")
SET(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -g -O0 -ftree-ter")
SET(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -g -O3 -Wstrict-aliasing=0")
SET(CMAKE_EXE_LINKER_FLAGS "")

SET(BUILD_FOR_XEON_PHI OFF)
MARK_AS_ADVANCED(BUILD_FOR_XEON_PHI)

SET(EXT "")

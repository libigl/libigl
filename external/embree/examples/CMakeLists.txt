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

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/common)

ADD_SUBDIRECTORY(common)
ADD_SUBDIRECTORY(renderer)
ADD_SUBDIRECTORY(tutorials)

SET(BUILD_TUTORIALS_SSE off CACHE BOOL "Enable to build Embree tutorials for SSE")

IF (BUILD_TUTORIALS_SSE)
  SET(TARGET_SPMD "sse")
  ADD_SUBDIRECTORY(tutorial00 tutorial00_sse)
  ADD_SUBDIRECTORY(tutorial01 tutorial01_sse)
  ADD_SUBDIRECTORY(tutorial02 tutorial02_sse)
  ADD_SUBDIRECTORY(tutorial03 tutorial03_sse)
ENDIF()

SET(BUILD_TUTORIALS_AVX off CACHE BOOL "Enable to build Embree tutorials for AVX")

IF (BUILD_TUTORIALS_AVX)
  SET(TARGET_SPMD "avx")
  ADD_SUBDIRECTORY(tutorial00 tutorial00_avx)
  ADD_SUBDIRECTORY(tutorial01 tutorial01_avx)
  ADD_SUBDIRECTORY(tutorial02 tutorial02_avx)
  ADD_SUBDIRECTORY(tutorial03 tutorial03_avx)
ENDIF()

SET(BUILD_TUTORIALS_KNC off CACHE BOOL "Enable to build Embree tutorials for KNC")

IF (BUILD_TUTORIALS_KNC)

  SET(TARGET_SPMD "knc")
  ADD_SUBDIRECTORY(tutorial_device)
  ADD_SUBDIRECTORY(tutorial_host) 

  SET (__HOST__ 1)
  ADD_SUBDIRECTORY(tutorial00 tutorial00_host_knc)
  ADD_SUBDIRECTORY(tutorial01 tutorial01_host_knc)
  ADD_SUBDIRECTORY(tutorial02 tutorial02_host_knc)
  ADD_SUBDIRECTORY(tutorial03 tutorial03_host_knc)

  SET (__HOST__ 0)
  ADD_SUBDIRECTORY(tutorial00 tutorial00_device_knc)
  ADD_SUBDIRECTORY(tutorial01 tutorial01_device_knc)
  ADD_SUBDIRECTORY(tutorial02 tutorial02_device_knc)
  ADD_SUBDIRECTORY(tutorial03 tutorial03_device_knc)
    
ENDIF()

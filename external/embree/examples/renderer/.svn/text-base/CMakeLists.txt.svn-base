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

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

ADD_SUBDIRECTORY(renderer)
ADD_SUBDIRECTORY(device)

SET(BUILD_NETWORK_DEVICE on CACHE BOOL "Enable to build network device")

IF (BUILD_NETWORK_DEVICE)
  SET(TARGET_CPU "xeon")
  ADD_SUBDIRECTORY(device_network)
ENDIF (BUILD_NETWORK_DEVICE)

SET(BUILD_SINGLERAY_DEVICE on CACHE BOOL "Enable to build single ray device")

IF (BUILD_SINGLERAY_DEVICE)
  SET(TARGET_CPU "xeon")
  ADD_SUBDIRECTORY(device_singleray)
ENDIF (BUILD_SINGLERAY_DEVICE)

SET(BUILD_SINGLERAY_DEVICE_KNC off CACHE BOOL "Enable to build single ray device for KNC")

IF (BUILD_SINGLERAY_DEVICE_KNC)
  SET(TARGET_CPU "xeon_phi")
  ADD_SUBDIRECTORY(device_singleray device_singleray_knc)
ENDIF (BUILD_SINGLERAY_DEVICE_KNC)

SET(BUILD_ISPC_DEVICE_SSE off CACHE BOOL "Enable to build ISPC SSE device")

IF (BUILD_ISPC_DEVICE_SSE)
  SET(TARGET_SPMD "sse")
  ADD_SUBDIRECTORY(device_ispc device_ispc_sse)
ENDIF (BUILD_ISPC_DEVICE_SSE)

SET(BUILD_ISPC_DEVICE_AVX off CACHE BOOL "Enable to build ISPC AVX device")

IF (BUILD_ISPC_DEVICE_AVX)
  SET(TARGET_SPMD "avx")
  ADD_SUBDIRECTORY(device_ispc device_ispc_avx)
ENDIF (BUILD_ISPC_DEVICE_AVX)

SET(BUILD_ISPC_DEVICE_KNC off CACHE BOOL "Enable to build ISPC device for KNC")

IF (BUILD_ISPC_DEVICE_KNC)
  SET(TARGET_SPMD "knc")
  ADD_SUBDIRECTORY(device_ispc device_ispc_knc)
ENDIF (BUILD_ISPC_DEVICE_KNC)

IF (BUILD_SINGLERAY_DEVICE_KNC OR BUILD_ISPC_DEVICE_KNC)
  SET(TARGET_CPU "xeon")
  ADD_SUBDIRECTORY(device_coi)
  SET(TARGET_CPU "xeon_phi")
  ADD_SUBDIRECTORY(device_coi device_coi_knc)
ENDIF ()

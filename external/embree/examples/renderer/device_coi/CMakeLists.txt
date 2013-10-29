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

CONFIGURE_CPU ()

FIND_PACKAGE(COI REQUIRED)
INCLUDE_DIRECTORIES(${COI_INCLUDE_PATHS})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF (__XEON_PHI__)

  ADD_LIBRARY(coi_server${EXT} STATIC coi_server.cpp)
  TARGET_LINK_LIBRARIES(coi_server${EXT} lexers${EXT} ${COI_DEV_LIBRARIES})

ELSE ()

  ADD_DEFINITIONS(-DNETWORK_NO_CREATE)
  ADD_LIBRARY(device_coi SHARED coi_device.cpp)
  TARGET_LINK_LIBRARIES(device_coi sys image ${COI_HOST_LIBRARIES})  

ENDIF ()

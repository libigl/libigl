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

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET (SOURCES 
    api/singleray_device.cpp
    lights/hdrilight.cpp   
    shapes/trianglemesh_normals.cpp   
    shapes/trianglemesh_full.cpp       
    samplers/sampler.cpp
    samplers/distribution1d.cpp
    samplers/distribution2d.cpp
    integrators/pathtraceintegrator.cpp
    filters/filter.cpp
    renderers/debugrenderer.cpp
    renderers/integratorrenderer.cpp
    renderers/progress.cpp
    )

IF (__XEON_PHI__)
  ADD_EXECUTABLE(device_singleray${EXT} ${SOURCES})
  TARGET_LINK_LIBRARIES(device_singleray${EXT} sys${EXT} embree${EXT} coi_server${EXT})
ELSE ()
  ADD_LIBRARY(device_singleray${EXT} SHARED ${SOURCES})
  TARGET_LINK_LIBRARIES(device_singleray${EXT} sys${EXT} embree${EXT} image${EXT})
ENDIF ()

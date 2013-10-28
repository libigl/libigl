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

CONFIGURE_SPMD()

SET (SOURCES 
  api/ispc_device.cpp
  shapes/trianglemesh.cpp
  shapes/sphere.cpp)

SET (ISPC_SOURCES
  api/ref.ispc
  scene/scene.ispc
  scene/instance.ispc
  cameras/pinholecamera.ispc
  cameras/depthoffieldcamera.ispc
  materials/material.ispc
  materials/matte.ispc
  materials/matte_textured.ispc
  materials/obj.ispc
  materials/dielectric.ispc
  materials/thindielectric.ispc
  materials/mirror.ispc
  materials/metal.ispc
  materials/velvet.ispc
  materials/metallicpaint.ispc
  materials/plastic.ispc
  lights/light.ispc
  lights/ambientlight.ispc
  lights/pointlight.ispc
  lights/directionallight.ispc
  lights/distantlight.ispc
  lights/spotlight.ispc
  lights/trianglelight.ispc
  lights/hdrilight.ispc
  shapes/shape.ispc
  shapes/trianglemesh.ispc
  tonemappers/defaulttonemapper.ispc
  textures/nearestneighbor.ispc
  textures/image3c.ispc
  textures/image3ca.ispc
  textures/image3f.ispc
  textures/image3fa.ispc
  renderers/renderer.ispc
  renderers/debugrenderer.ispc
  renderers/pathtracer.ispc
  framebuffers/swapchain.ispc
  framebuffers/accubuffer.ispc
  framebuffers/framebuffer.ispc
  framebuffers/framebuffer_rgb_float32.ispc
  framebuffers/framebuffer_rgba8.ispc
  samplers/distribution2d.ispc
  )

IF (__XEON__)
  SET (ISPC_SOURCES ${ISPC_SOURCES} 
    framebuffers/framebuffer_rgb8.ispc)
ENDIF ()

INCLUDE (build_ispc)

IF (__XEON_PHI__)
  ADD_EXECUTABLE(device_ispc${EXT} ${SOURCES})
  TARGET_LINK_LIBRARIES(device_ispc${EXT} sys${EXT} embree_ispc${EXT} coi_server${EXT})
ELSE (__XEON_PHI__)
  ADD_LIBRARY(device_ispc${EXT} SHARED ${SOURCES})
  TARGET_LINK_LIBRARIES(device_ispc${EXT} sys image embree_ispc${EXT})
ENDIF (__XEON_PHI__)

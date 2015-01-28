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

IF (__XEON__ AND BUILD_TUTORIALS)
  ADD_EXECUTABLE(tutorial05 tutorial05.cpp tutorial05_device.cpp)
  TARGET_LINK_LIBRARIES(tutorial05 embree tutorial transport tutorial_device)
ENDIF()

IF (__XEON__ AND BUILD_TUTORIALS_ISPC)
  ADD_ISPC_EXECUTABLE(tutorial05_ispc tutorial05.cpp tutorial05_device.ispc)
  TARGET_LINK_LIBRARIES(tutorial05_ispc embree tutorial transport tutorial_device_ispc)
ENDIF()

IF (NOT __XEON__ AND __HOST__)
  ADD_EXECUTABLE(tutorial05_xeonphi tutorial05.cpp)
  TARGET_LINK_LIBRARIES(tutorial05_xeonphi tutorial transport_host)
ENDIF()

IF (NOT __XEON__ AND NOT __HOST__)
  INCLUDE (icc_xeonphi)
  ADD_ISPC_EXECUTABLE(tutorial05_xeonphi_device tutorial05_device.ispc)
  TARGET_LINK_LIBRARIES(tutorial05_xeonphi_device embree_xeonphi.so transport_device tutorial_xeonphi_device_ispc)
  ADD_DEPENDENCIES(tutorial05_xeonphi_device embree_xeonphi)
ENDIF()

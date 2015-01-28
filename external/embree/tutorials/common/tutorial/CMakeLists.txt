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

IF (__XEON__)

FIND_PACKAGE(GLUT REQUIRED)
FIND_PACKAGE(OpenGL REQUIRED)

ADD_LIBRARY(tutorial STATIC 
    glutdisplay.cpp
    obj_loader.cpp)
TARGET_LINK_LIBRARIES(tutorial sys lexers ${OPENGL_LIBRARY} ${GLUT_LIBRARY})

IF (BUILD_TUTORIALS)
  ADD_LIBRARY(tutorial_device STATIC tutorial_device.cpp)
ENDIF()

IF (BUILD_TUTORIALS_ISPC)
  ADD_ISPC_LIBRARY(tutorial_device_ispc STATIC tutorial_device.ispc)
  SET_TARGET_PROPERTIES(tutorial_device_ispc PROPERTIES LINKER_LANGUAGE C)
ENDIF()

ELSE ()

IF (BUILD_TUTORIALS_ISPC)
  INCLUDE (icc_xeonphi)
  ADD_ISPC_LIBRARY(tutorial_xeonphi_device_ispc STATIC tutorial_device.ispc)
  SET_TARGET_PROPERTIES(tutorial_xeonphi_device_ispc PROPERTIES LINKER_LANGUAGE C)
ENDIF()

ENDIF ()

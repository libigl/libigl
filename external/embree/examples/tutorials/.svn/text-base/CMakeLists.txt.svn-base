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

FIND_PACKAGE(GLUT REQUIRED)
FIND_PACKAGE(OpenGL REQUIRED)

ADD_LIBRARY(tutorials STATIC 
    glutdisplay.cpp
    obj_loader.cpp)
TARGET_LINK_LIBRARIES(tutorials sys lexers ${OPENGL_LIBRARY} ${GLUT_LIBRARY})

ADD_LIBRARY(tutorial_host_device STATIC tutorials_host_device.cpp)
TARGET_LINK_LIBRARIES(tutorial_host_device tutorials ${OPENGL_LIBRARY} ${GLUT_LIBRARY})

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

SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_MODULE_PATH})

FIND_PACKAGE(OpenEXR)
IF(OPENEXR_FOUND)
  ADD_DEFINITIONS(-DUSE_OPENEXR)
  INCLUDE_DIRECTORIES(${OPENEXR_INCLUDE_PATHS})
ENDIF(OPENEXR_FOUND)

FIND_PACKAGE(ImageMagick)
IF(IMAGEMAGICK_FOUND)
  ADD_DEFINITIONS(-DUSE_IMAGEMAGICK)
  INCLUDE_DIRECTORIES(${IMAGEMAGICK_INCLUDE_PATHS})
ENDIF(IMAGEMAGICK_FOUND)

ADD_LIBRARY(image STATIC
  image.cpp
  ppm.cpp
  pfm.cpp
  tga.cpp
  exr.cpp
  magick.cpp
)
TARGET_LINK_LIBRARIES(image sys ${OPENEXR_LIBRARIES} ${IMAGEMAGICK_LIBRARIES})

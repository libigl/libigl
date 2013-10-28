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

SET(LIBRARY_PATHS
  /usr/lib /usr/local/lib /opt/local/lib)

FIND_PATH(IMAGEMAGICK_INCLUDE_PATH Magick++.h
  /usr/include       /usr/include/ImageMagick
  /usr/local/include /usr/local/include/ImageMagick
  /opt/local/include /opt/local/include/ImageMagick
)

FIND_LIBRARY(IMAGEMAGICK_CXX_LIBRARY NAMES Magick++          PATHS ${LIBRARY_PATHS} )
FIND_LIBRARY(IMAGEMAGICK_C_LIBRARY   NAMES Magick MagickCore PATHS ${LIBRARY_PATHS} )

IF (IMAGEMAGICK_INCLUDE_PATH AND IMAGEMAGICK_CXX_LIBRARY AND IMAGEMAGICK_C_LIBRARY)
  SET(IMAGEMAGICK_FOUND TRUE)
  SET(IMAGEMAGICK_INCLUDE_PATHS ${IMAGEMAGICK_INCLUDE_PATH} CACHE STRING "ImageMagick include paths")
  SET(IMAGEMAGICK_LIBRARIES ${IMAGEMAGICK_CXX_LIBRARY} ${IMAGEMAGICK_C_LIBRARY} CACHE STRING "ImageMagick Libraries")
ENDIF (IMAGEMAGICK_INCLUDE_PATH AND IMAGEMAGICK_CXX_LIBRARY AND IMAGEMAGICK_C_LIBRARY)

MARK_AS_ADVANCED(
  IMAGEMAGICK_FOUND
  IMAGEMAGICK_INCLUDE_PATHS
  IMAGEMAGICK_LIBRARIES
  IMAGEMAGICK_CXX_LIBRARY
  IMAGEMAGICK_C_LIBRARY
  IMAGEMAGICK_INCLUDE_PATH
)

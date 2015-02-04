## ======================================================================== ##
## Copyright 2009-2014 Intel Corporation                                    ##
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
INCLUDE(CMakeDependentOption)

SET(ADDITIONAL_SOURCES)

FIND_PACKAGE(ImageMagick COMPONENTS Magick++)
MARK_AS_ADVANCED(
  ImageMagick_EXECUTABLE_DIR
  ImageMagick_Magick++_INCLUDE_DIR
  ImageMagick_Magick++_LIBRARY
)
CMAKE_DEPENDENT_OPTION(USE_IMAGE_MAGICK "Enables BMP, GIF, PNG, TGA, TIFF image codecs." OFF "IMAGEMAGICK_FOUND;BUILD_TUTORIALS OR BUILD_TUTORIALS_ISPC" OFF)
IF (USE_IMAGE_MAGICK)
  ADD_DEFINITIONS(-DUSE_IMAGEMAGICK)
  INCLUDE_DIRECTORIES(${ImageMagick_Magick++_INCLUDE_DIR})
  SET(ADDITIONAL_LIBRARIES ${ADDITIONAL_LIBRARIES} ${ImageMagick_Magick++_LIBRARY})
  SET(ADDITIONAL_SOURCES ${ADDITIONAL_SOURCES} magick.cpp)
ENDIF (USE_IMAGE_MAGICK)

FIND_PACKAGE(JPEG)
CMAKE_DEPENDENT_OPTION(USE_LIBJPEG "Enables JPEG image codec." ON "JPEG_FOUND;BUILD_TUTORIALS OR BUILD_TUTORIALS_ISPC" OFF)
IF (USE_LIBJPEG)
  ADD_DEFINITIONS(-DUSE_LIBJPEG)
  INCLUDE_DIRECTORIES(${JPEG_INCLUDE_DIR})
  SET(ADDITIONAL_LIBRARIES ${ADDITIONAL_LIBRARIES} ${JPEG_LIBRARIES})
  SET(ADDITIONAL_SOURCES ${ADDITIONAL_SOURCES} jpeg.cpp)
ENDIF (USE_LIBJPEG)

FIND_PACKAGE(OpenEXR)
CMAKE_DEPENDENT_OPTION(USE_OPENEXR "Enables OpenEXR image codec." ON "OPENEXR_FOUND;BUILD_TUTORIALS OR BUILD_TUTORIALS_ISPC" OFF)
IF (USE_OPENEXR)
  ADD_DEFINITIONS(-DUSE_OPENEXR)
  INCLUDE_DIRECTORIES(${OPENEXR_INCLUDE_DIRS})
  SET(ADDITIONAL_LIBRARIES ${ADDITIONAL_LIBRARIES} ${OPENEXR_LIBRARIES})
  SET(ADDITIONAL_SOURCES ${ADDITIONAL_SOURCES} exr.cpp)
ENDIF (USE_OPENEXR)

IF (BUILD_TUTORIALS OR BUILD_TUTORIALS_ISPC)
  ADD_LIBRARY(image STATIC
    image.cpp
    pfm.cpp
    ppm.cpp
    tga.cpp
    ${ADDITIONAL_SOURCES}
  )
  TARGET_LINK_LIBRARIES(image sys ${ADDITIONAL_LIBRARIES})
  SET_PROPERTY(TARGET image PROPERTY FOLDER tutorials/common)
ENDIF ()

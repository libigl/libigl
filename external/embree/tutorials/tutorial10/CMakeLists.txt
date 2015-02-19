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

INCLUDE(CMakeDependentOption)

FIND_PACKAGE(OpenSubdiv REQUIRED)
CMAKE_DEPENDENT_OPTION(USE_OPENSUBDIV "Uses OpenSubdiv for verification in tutorial10." ON "OPENSUBDIV_FOUND;BUILD_TUTORIALS" OFF)
SET(ADDITIONAL_LIBRARIES)
IF (USE_OPENSUBDIV)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${OPENSUBDIV_INCLUDE_DIRS})
  ADD_DEFINITIONS(-D__USE_OPENSUBDIV__)
  SET(ADDITIONAL_LIBRARIES ${ADDITIONAL_LIBRARIES} ${OPENSUBDIV_LIBRARIES})
ENDIF()

INCLUDE(tutorial)
ADD_TUTORIAL(tutorial10)

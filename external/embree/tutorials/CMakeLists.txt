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

INCLUDE_DIRECTORIES     (${CMAKE_CURRENT_SOURCE_DIR}/common)
INCLUDE_DIRECTORIES_ISPC(${CMAKE_CURRENT_SOURCE_DIR}/common)

ADD_SUBDIRECTORY(common)

SET(BUILD_TUTORIALS      ON CACHE BOOL "Enable to build Embree tutorials")
SET(BUILD_TUTORIALS_ISPC ON CACHE BOOL "Enable to build Embree ISPC tutorials")

IF (BUILD_TUTORIALS OR BUILD_TUTORIALS_ISPC)
  SET (__XEON__ 1)
  ADD_SUBDIRECTORY(tutorial00)
  ADD_SUBDIRECTORY(tutorial01)
  ADD_SUBDIRECTORY(tutorial02)
  ADD_SUBDIRECTORY(tutorial03)
  ADD_SUBDIRECTORY(tutorial04)
  ADD_SUBDIRECTORY(tutorial05)
  ADD_SUBDIRECTORY(tutorial06)
ENDIF()

IF (TARGET_XEON_PHI)
  SET (__XEON__ 0)
  SET (__HOST__ 1)
  ADD_SUBDIRECTORY(tutorial00 tutorial00_host)
  ADD_SUBDIRECTORY(tutorial01 tutorial01_host)
  ADD_SUBDIRECTORY(tutorial02 tutorial02_host)
  ADD_SUBDIRECTORY(tutorial03 tutorial03_host)
  ADD_SUBDIRECTORY(tutorial04 tutorial04_host)
  ADD_SUBDIRECTORY(tutorial05 tutorial05_host)
  ADD_SUBDIRECTORY(tutorial06 tutorial06_host)
  SET (__HOST__ 0)
  ADD_SUBDIRECTORY(tutorial00 tutorial00_device)
  ADD_SUBDIRECTORY(tutorial01 tutorial01_device)
  ADD_SUBDIRECTORY(tutorial02 tutorial02_device)
  ADD_SUBDIRECTORY(tutorial03 tutorial03_device)
  ADD_SUBDIRECTORY(tutorial04 tutorial04_device)
  ADD_SUBDIRECTORY(tutorial05 tutorial05_device)
  ADD_SUBDIRECTORY(tutorial06 tutorial06_device)
ENDIF()

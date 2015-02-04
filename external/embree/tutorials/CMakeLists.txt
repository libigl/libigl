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

INCLUDE_DIRECTORIES     (${CMAKE_CURRENT_SOURCE_DIR}/common)
INCLUDE_DIRECTORIES_ISPC(${CMAKE_CURRENT_SOURCE_DIR}/common)

OPTION(BUILD_TUTORIALS      "Enable to build Embree tutorials" ON)
OPTION(BUILD_TUTORIALS_ISPC "Enable to build Embree ISPC tutorials" ON)

ADD_SUBDIRECTORY(common)

SET (__XEON__ 1)
ADD_SUBDIRECTORY(tutorial00)
ADD_SUBDIRECTORY(tutorial01)
ADD_SUBDIRECTORY(tutorial02)
ADD_SUBDIRECTORY(tutorial03)
ADD_SUBDIRECTORY(tutorial04)
ADD_SUBDIRECTORY(tutorial05)
ADD_SUBDIRECTORY(tutorial06)
ADD_SUBDIRECTORY(tutorial07)
ADD_SUBDIRECTORY(tutorial08)
ADD_SUBDIRECTORY(tutorial09)
ADD_SUBDIRECTORY(tutorial10)

IF (XEON_PHI_ISA)
  SET (__XEON__ 0)
  SET (__HOST__ 1)
  ADD_SUBDIRECTORY(tutorial00 tutorial00_host)
  ADD_SUBDIRECTORY(tutorial01 tutorial01_host)
  ADD_SUBDIRECTORY(tutorial02 tutorial02_host)
  ADD_SUBDIRECTORY(tutorial03 tutorial03_host)
  ADD_SUBDIRECTORY(tutorial04 tutorial04_host)
  ADD_SUBDIRECTORY(tutorial05 tutorial05_host)
  ADD_SUBDIRECTORY(tutorial06 tutorial06_host)
  ADD_SUBDIRECTORY(tutorial07 tutorial07_host)
  ADD_SUBDIRECTORY(tutorial08 tutorial08_host)
  ADD_SUBDIRECTORY(tutorial09 tutorial09_host)
  ADD_SUBDIRECTORY(tutorial10 tutorial10_host)
  SET (__HOST__ 0)
  ADD_SUBDIRECTORY(tutorial00 tutorial00_device)
  ADD_SUBDIRECTORY(tutorial01 tutorial01_device)
  ADD_SUBDIRECTORY(tutorial02 tutorial02_device)
  ADD_SUBDIRECTORY(tutorial03 tutorial03_device)
  ADD_SUBDIRECTORY(tutorial04 tutorial04_device)
  ADD_SUBDIRECTORY(tutorial05 tutorial05_device)
  ADD_SUBDIRECTORY(tutorial06 tutorial06_device)
  ADD_SUBDIRECTORY(tutorial07 tutorial07_device)
  ADD_SUBDIRECTORY(tutorial08 tutorial08_device)
  ADD_SUBDIRECTORY(tutorial09 tutorial09_device)
  ADD_SUBDIRECTORY(tutorial10 tutorial10_device)
ENDIF()

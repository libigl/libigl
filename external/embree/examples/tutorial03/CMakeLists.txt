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

SET (HOST_SOURCES tutorial03.cpp)
SET (ISPC_SOURCES tutorial03.ispc)

IF (__HOST__)

  ADD_EXECUTABLE(tutorial03_knc ${HOST_SOURCES})
  TARGET_LINK_LIBRARIES(tutorial03_knc tutorials_host)

ELSE ()

  CONFIGURE_SPMD()
  INCLUDE (build_ispc)

  IF (__XEON__)
    ADD_EXECUTABLE(tutorial03${EXT} ${SOURCES} ${HOST_SOURCES})
    TARGET_LINK_LIBRARIES(tutorial03${EXT} tutorial_host_device embree_ispc${EXT})
  ELSE ()
    ADD_EXECUTABLE(tutorial03_device${EXT} ${SOURCES})
    TARGET_LINK_LIBRARIES(tutorial03_device${EXT} tutorials_device${EXT} embree_ispc${EXT})
  ENDIF ()

ENDIF ()


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

SET (HOST_SOURCES tutorial02.cpp)
SET (ISPC_SOURCES tutorial02.ispc)

IF (__HOST__)

  ADD_EXECUTABLE(tutorial02_knc ${HOST_SOURCES})
  TARGET_LINK_LIBRARIES(tutorial02_knc tutorials_host)

ELSE ()

  CONFIGURE_SPMD()
  INCLUDE (build_ispc)

  IF (__XEON__)
    ADD_EXECUTABLE(tutorial02${EXT} ${SOURCES} ${HOST_SOURCES})
    TARGET_LINK_LIBRARIES(tutorial02${EXT} tutorial_host_device embree_ispc${EXT})
  ELSE ()
    ADD_EXECUTABLE(tutorial02_device${EXT} ${SOURCES})
    TARGET_LINK_LIBRARIES(tutorial02_device${EXT} tutorials_device${EXT} embree_ispc${EXT})
  ENDIF ()

ENDIF ()


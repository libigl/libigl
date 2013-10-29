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

CONFIGURE_CPU ()

ADD_LIBRARY(sys${EXT} STATIC
  platform.cpp
  sysinfo.cpp
  filename.cpp
  library.cpp
  thread.cpp
  network.cpp
  taskscheduler.cpp
  taskscheduler_sys.cpp
  taskscheduler_mic.cpp
  sync/mutex.cpp
  sync/condition.cpp
  stl/string.cpp
)

IF (__XEON_PHI__)
  TARGET_LINK_LIBRARIES(sys${EXT} simd${EXT} pthread dl)
ELSE ()
  TARGET_LINK_LIBRARIES(sys${EXT} pthread dl)
ENDIF() 

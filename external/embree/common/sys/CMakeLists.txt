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

SET(CMAKE_THREAD_PREFER_PTHREAD TRUE)
FIND_PACKAGE(Threads REQUIRED)

IF (__XEON__)

IF (TARGET_AVX)
  ADD_DEFINITIONS(-D__TARGET_AVX__)
ENDIF()

IF (TARGET_AVX2)
  ADD_DEFINITIONS(-D__TARGET_AVX2__)
ENDIF()

ADD_LIBRARY(sys STATIC
  platform.cpp
  sysinfo.cpp
  filename.cpp
  library.cpp
  thread.cpp
  network.cpp
  tasklogger.cpp
  taskscheduler.cpp
  taskscheduler_sys.cpp
  sync/mutex.cpp
  sync/condition.cpp
  sync/barrier.cpp
  stl/string.cpp
)

TARGET_LINK_LIBRARIES(sys ${CMAKE_THREAD_LIBS_INIT} ${CMAKE_DL_LIBS})
SET_PROPERTY(TARGET sys PROPERTY FOLDER kernels)

ELSE ()
 
INCLUDE(icc_xeonphi)

ADD_LIBRARY(sys_xeonphi STATIC
  platform.cpp
  sysinfo.cpp
  filename.cpp
  library.cpp
  thread.cpp
  network.cpp
  tasklogger.cpp
  taskscheduler.cpp
  taskscheduler_sys.cpp
  taskscheduler_mic.cpp
  sync/mutex.cpp
  sync/condition.cpp
  sync/barrier.cpp
  stl/string.cpp
)

TARGET_LINK_LIBRARIES(sys_xeonphi ${CMAKE_THREAD_LIBS_INIT} ${CMAKE_DL_LIBS})

ENDIF ()

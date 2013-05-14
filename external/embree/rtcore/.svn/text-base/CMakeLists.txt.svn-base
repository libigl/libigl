## ======================================================================== ##
## Copyright 2009-2012 Intel Corporation                                    ##
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

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(USE_STAT_COUNTERS false CACHE BOOL "Set to 1 to activate statistics counters")

IF (USE_STAT_COUNTERS)
ADD_DEFINITIONS(-D__USE_STAT_COUNTERS__)
ENDIF (USE_STAT_COUNTERS)

ADD_LIBRARY(rtcore STATIC

  common/accel.cpp
  common/alloc.cpp
  common/stat.cpp 
  common/primrefgen.cpp 
  common/heuristic_binning.cpp 
  common/heuristic_spatial.cpp 
  common/splitter.cpp 
  common/splitter_parallel.cpp 
  common/splitter_fallback.cpp 

  bvh2/bvh2.cpp   
  bvh2/bvh2_intersector.cpp   
  bvh2/bvh2_builder.cpp   

  bvh4/bvh4.cpp   
  bvh4/bvh4_intersector.cpp   
  bvh4/bvh4_builder.cpp   

  bvh4mb/bvh4mb.cpp   
  bvh4mb/bvh4mb_builder.cpp   
  bvh4mb/bvh4mb_intersector.cpp   
)

TARGET_LINK_LIBRARIES(rtcore sys)

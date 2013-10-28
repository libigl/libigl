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

CONFIGURE_CPU()

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(USE_STAT_COUNTERS false CACHE BOOL "Set to 1 to activate statistics counters")

IF (USE_STAT_COUNTERS)
  ADD_DEFINITIONS(-D__USE_STAT_COUNTERS__)
ENDIF (USE_STAT_COUNTERS)

ADD_SUBDIRECTORY(sys)
ADD_SUBDIRECTORY(simd)

SET(SOURCES
  embree.cpp
  common/alloc.cpp
  common/stat.cpp 
  common/accel.cpp
  common/registry_accel.cpp
  common/registry_builder.cpp
  common/registry_intersector.cpp

  geometry/triangle_mesh.cpp
  geometry/virtual_scene.cpp

  builders/primrefgen.cpp 
  builders/heuristic_binning.cpp 
  builders/heuristic_spatial.cpp 
  builders/splitter.cpp 
  builders/splitter_parallel.cpp 
  builders/splitter_fallback.cpp 

  bvh2/bvh2.cpp 
)

IF (__XEON__)

  SET(SOURCES ${SOURCES} 

  
  bvh2/bvh2_intersector1.cpp   
  bvh2/bvh2_intersector4.cpp   
  bvh2/bvh2_intersector8.cpp   
  
  bvh4/bvh4.cpp   
  bvh4/bvh4_intersector1.cpp   
  bvh4/bvh4_intersector4_single.cpp   
  bvh4/bvh4_intersector4_chunk.cpp   
  bvh4/bvh4_intersector4_hybrid.cpp   
  bvh4/bvh4_intersector8_single.cpp   
  bvh4/bvh4_intersector8_chunk.cpp   
  bvh4/bvh4_intersector8_hybrid.cpp   

  bvh4i/bvh4i.cpp   
  bvh4i/bvh4i_intersector1.cpp   
                           
  bvh4mb/bvh4mb.cpp
  bvh4mb/bvh4mb_builder.cpp   
  bvh4mb/bvh4mb_intersector1.cpp   
  bvh4mb/bvh4mb_intersector4.cpp   
  bvh4mb/bvh4mb_intersector8.cpp
  )

  IF (__AVX__)
    SET(SOURCES ${SOURCES} 
    bvh4/bvh4_intersector1_avx.cpp
    bvh8/bvh8.cpp   
    bvh8/bvh8_intersector1.cpp 
    )
  ENDIF ()

ENDIF ()

IF (__XEON_PHI__)

  SET(SOURCES ${SOURCES}
  bvh4aos/bvh4aos.cpp   

  bvh4aos/bvh4aos_intersector1.cpp   
  bvh4aos/bvh4aos_intersector16_chunk.cpp   

  bvh4aos/bvh4aos_triangle1_intersector1.cpp   
  bvh4aos/bvh4aos_triangle1_intersector16_single.cpp   
  bvh4aos/bvh4aos_triangle1_intersector16_chunk.cpp   
  bvh4aos/bvh4aos_triangle1_intersector16_hybrid.cpp   

  bvh4aos/bvh4aos_triangle1_intersector1_ref.cpp
  bvh4aos/bvh4aos_triangle1_intersector16_chunk_ref.cpp   
  bvh4aos/bvh4aos_triangle1_intersector16_hybrid_ref.cpp   

  bvh4aos/builder/bvh4aos_builder.cpp   
  bvh4aos/builder/bvh4aos_task_scheduler.cpp   
  bvh4aos/builder/bvh4aos_builder_util.cpp 
  bvh4aos/builder/bvh4aos_builder_sah.cpp 
  bvh4aos/builder/bvh4aos_builder_morton_code.cpp   

  )
ENDIF()

ADD_LIBRARY(embree${EXT} STATIC ${SOURCES})
TARGET_LINK_LIBRARIES(embree${EXT} sys${EXT})

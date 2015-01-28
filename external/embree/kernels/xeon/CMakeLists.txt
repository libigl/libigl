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

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF (TARGET_SSE41)
 ADD_DEFINITIONS(-D__TARGET_SSE41__)
ENDIF()

IF (TARGET_AVX)
 ADD_DEFINITIONS(-D__TARGET_AVX__)
ENDIF()

IF (TARGET_AVX2)
 ADD_DEFINITIONS(-D__TARGET_AVX2__)
ENDIF()

ADD_ISPC_LIBRARY(embree SHARED   

  ../common/stat.cpp 
  ../common/alloc.cpp 
  ../common/tasksys.cpp 
  ../common/acceln.cpp
  ../common/rtcore.cpp 
  ../common/rtcore_ispc.cpp 
  ../common/rtcore_ispc.ispc 
  ../common/buffer.cpp
  ../common/scene.cpp
  ../common/geometry.cpp
  ../common/scene_user_geometry.cpp
  ../common/scene_triangle_mesh.cpp
  ../common/scene_quadratic_bezier_curves.cpp
  
  builders/heuristic_binning.cpp
  builders/heuristic_spatial.cpp
  builders/splitter.cpp
  builders/splitter_fallback.cpp
  builders/splitter_parallel.cpp
  builders/primrefgen.cpp
  
  geometry/triangle1.cpp
  geometry/triangle4.cpp
  geometry/triangle1v.cpp
  geometry/triangle4v.cpp
  geometry/triangle4i.cpp
  geometry/ispc_wrapper_sse.cpp
  geometry/instance_intersector1.cpp
  geometry/instance_intersector4.cpp
  
  bvh4/bvh4.cpp
  bvh4/bvh4_rotate.cpp
  bvh4/bvh4_refit.cpp
  bvh4/bvh4_builder.cpp
  bvh4/bvh4_builder_fast.cpp
  bvh4/bvh4_builder_morton.cpp
  bvh4/bvh4_builder_binner.cpp
  bvh4/bvh4_builder_toplevel.cpp
  bvh4/bvh4_intersector1.cpp   
  bvh4/bvh4_intersector4_chunk.cpp
  bvh4/bvh4_intersector4_hybrid.cpp
  bvh4/bvh4_statistics.cpp
  bvh4/virtual_accel.cpp
  bvh4/twolevel_accel.cpp

  bvh4i/bvh4i.cpp
  bvh4i/bvh4i_statistics.cpp
  bvh4i/bvh4i_rotate.cpp
  bvh4i/bvh4i_builder.cpp
  bvh4i/bvh4i_builder_binner.cpp
  bvh4i/bvh4i_intersector1.cpp   
  bvh4i/bvh4i_intersector4_chunk.cpp   
  bvh4i/bvh4i_intersector1_scalar.cpp   
 
  bvh4mb/bvh4mb.cpp
  bvh4mb/bvh4mb_builder.cpp
  bvh4mb/bvh4mb_intersector1.cpp   
  bvh4mb/bvh4mb_intersector4.cpp
)

IF (TARGET_SSE41) 
  ADD_LIBRARY(embree_sse41 STATIC
    bvh4/bvh4_builder_binner.cpp
    bvh4/bvh4_builder_toplevel.cpp
    bvh4/bvh4_intersector1.cpp   
    bvh4/bvh4_intersector4_chunk.cpp
    bvh4/bvh4_intersector4_hybrid.cpp
    bvh4/bvh4_builder_morton.cpp
  )
  SET_TARGET_PROPERTIES(embree_sse41 PROPERTIES COMPILE_FLAGS "${FLAGS_SSE41}")
  SET(EMBREE_LIBRARIES ${EMBREE_LIBRARIES} embree_sse41)
ENDIF ()

IF (TARGET_AVX) 
  ADD_LIBRARY(embree_avx STATIC
   
   geometry/triangle8.cpp
   geometry/ispc_wrapper_avx.cpp

   geometry/instance_intersector1.cpp
   geometry/instance_intersector4.cpp
   geometry/instance_intersector8.cpp

   bvh4/bvh4_intersector1.cpp   
   bvh4/bvh4_intersector4_chunk.cpp
   bvh4/bvh4_intersector4_hybrid.cpp
   bvh4/bvh4_intersector8_chunk.cpp
   bvh4/bvh4_intersector8_hybrid.cpp

   bvh4i/bvh4i_intersector1.cpp   
   bvh4i/bvh4i_intersector1_scalar.cpp   
   bvh4i/bvh4i_intersector4_chunk.cpp   
   bvh4i/bvh4i_intersector8_chunk.cpp   
   bvh4i/bvh4i_intersector8_hybrid.cpp

   bvh4i/bvh4i_builder_fast.cpp
   bvh4i/bvh4i_builder_binner.cpp
   bvh4i/bvh4i_builder_morton.cpp
   bvh4i/bvh4i_builder_morton_enhanced.cpp

   bvh8i/bvh8i.cpp
   bvh8i/bvh8i_builder.cpp
   bvh8i/bvh8i_intersector1.cpp   
   bvh8i/bvh8i_intersector8_chunk.cpp
   bvh8i/bvh8i_intersector8_hybrid.cpp

   bvh4mb/bvh4mb_intersector1.cpp   
   bvh4mb/bvh4mb_intersector4.cpp
   bvh4mb/bvh4mb_intersector8.cpp
)

  SET_TARGET_PROPERTIES(embree_avx PROPERTIES COMPILE_FLAGS "${FLAGS_AVX}")
  SET(EMBREE_LIBRARIES ${EMBREE_LIBRARIES} embree_avx)
 ENDIF()

IF (TARGET_AVX2) 
  ADD_LIBRARY(embree_avx2 STATIC

    geometry/instance_intersector1.cpp
    geometry/instance_intersector4.cpp
    geometry/instance_intersector8.cpp

    bvh4/bvh4_intersector1.cpp
    bvh4/bvh4_intersector4_chunk.cpp
    bvh4/bvh4_intersector4_hybrid.cpp
    bvh4/bvh4_intersector8_chunk.cpp
    bvh4/bvh4_intersector8_hybrid.cpp
    bvh4i/bvh4i_intersector8_chunk_avx2.cpp  
    bvh4i/bvh4i_intersector8_hybrid.cpp

    bvh4i/bvh4i_intersector1.cpp   
    bvh4i/bvh4i_intersector1_scalar.cpp   
    bvh4i/bvh4i_intersector4_chunk.cpp   
    bvh4i/bvh4i_intersector8_chunk.cpp

    bvh4mb/bvh4mb_intersector1.cpp
    bvh4mb/bvh4mb_intersector4.cpp
    bvh4mb/bvh4mb_intersector8.cpp

    bvh8i/bvh8i_intersector1.cpp   
    bvh8i/bvh8i_intersector8_chunk.cpp
    bvh8i/bvh8i_intersector8_hybrid.cpp

    bvh4i/bvh4i_builder_fast.cpp
    bvh4i/bvh4i_builder_binner.cpp
)

  SET_TARGET_PROPERTIES(embree_avx2 PROPERTIES COMPILE_FLAGS "${FLAGS_AVX2}")
  SET(EMBREE_LIBRARIES ${EMBREE_LIBRARIES} embree_avx2)
ENDIF()

TARGET_LINK_LIBRARIES(embree ${EMBREE_LIBRARIES} sys simd)
SET_TARGET_PROPERTIES(embree PROPERTIES VERSION 2.2.0 SOVERSION 2)
INSTALL(TARGETS embree DESTINATION lib)

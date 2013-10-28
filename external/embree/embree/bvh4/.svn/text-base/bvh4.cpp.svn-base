// ======================================================================== //
// Copyright 2009-2013 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#include "bvh4.h"
#include "../geometry/triangles.h"
#include "../builders/heuristics.h"
#include "../builders/bvh_builder.h"
#include "../builders/bvh_statistics.h"

namespace embree
{
  /*! intersector registration functions */
  void BVH4Intersector1Register();
  void BVH4Intersector4SingleRegister();
  void BVH4Intersector4ChunkRegister();
  void BVH4Intersector4HybridRegister();
  void BVH4Intersector8SingleRegister();
  void BVH4Intersector8ChunkRegister();
  void BVH4Intersector8HybridRegister();
  void BVH4Intersector1AVXRegister();
  
  void BVH4Register () 
  {
    /* register acceleration structure */
    TriangleMesh::accels.add("bvh4",BVH4::create);
    VirtualScene   ::accels.add("bvh4",BVH4::create);
#if !defined (__MIC__)
    TriangleMesh::accels.setDefaultAccel("bvh4");
    VirtualScene   ::accels.setDefaultAccel("bvh4");
#endif
    TriangleMesh::accels.setDefaultTriangle("bvh4","triangle4");
    VirtualScene   ::accels.setDefaultTriangle("bvh4","virtual");
    TriangleMesh::builders.setDefaultBuilder ("bvh4","objectsplit");
    VirtualScene   ::builders.setDefaultBuilder ("bvh4","objectsplit");

    /* register triangle mesh builders */
    TriangleMesh::builders.add("bvh4","triangle1i","objectsplit" ,BVHBuilderT<BVH4, HeuristicBinning<0> >::create,1,inf);
    TriangleMesh::builders.add("bvh4","triangle1v","objectsplit" ,BVHBuilderT<BVH4, HeuristicBinning<0> >::create,1,inf);
    TriangleMesh::builders.add("bvh4","triangle4i","objectsplit" ,BVHBuilderT<BVH4, HeuristicBinning<2> >::create,1,inf);
    TriangleMesh::builders.add("bvh4","triangle4v","objectsplit" ,BVHBuilderT<BVH4, HeuristicBinning<2> >::create,1,inf);
    TriangleMesh::builders.add("bvh4","triangle1" ,"objectsplit", BVHBuilderT<BVH4, HeuristicBinning<0> >::create,1,inf);
    TriangleMesh::builders.add("bvh4","triangle4" ,"objectsplit", BVHBuilderT<BVH4, HeuristicBinning<2> >::create,1,inf);
#if defined (__AVX__)
    TriangleMesh::builders.add("bvh4","triangle8" ,"objectsplit", BVHBuilderT<BVH4, HeuristicBinning<3> >::create,1,inf);
#endif

    TriangleMesh::builders.add("bvh4","triangle1i","spatialsplit" ,BVHBuilderT<BVH4, HeuristicSpatial<0> >::create,1,inf);
    TriangleMesh::builders.add("bvh4","triangle1v","spatialsplit" ,BVHBuilderT<BVH4, HeuristicSpatial<0> >::create,1,inf);
    TriangleMesh::builders.add("bvh4","triangle4i","spatialsplit" ,BVHBuilderT<BVH4, HeuristicSpatial<2> >::create,1,inf);
    TriangleMesh::builders.add("bvh4","triangle4v","spatialsplit" ,BVHBuilderT<BVH4, HeuristicSpatial<2> >::create,1,inf);
    TriangleMesh::builders.add("bvh4","triangle1" ,"spatialsplit", BVHBuilderT<BVH4, HeuristicSpatial<0> >::create,1,inf);
    TriangleMesh::builders.add("bvh4","triangle4" ,"spatialsplit", BVHBuilderT<BVH4, HeuristicSpatial<2> >::create,1,inf);
#if defined (__AVX__)
    TriangleMesh::builders.add("bvh4","triangle8" ,"spatialsplit", BVHBuilderT<BVH4, HeuristicSpatial<3> >::create,1,inf);
#endif

    /* register virtual object builders */
    VirtualScene::builders.add("bvh4","virtual","objectsplit", BVHBuilderT<BVH4, HeuristicBinning<0> >::create,1,inf);

    /* register intersectors */
#if defined (__SSE__)
    BVH4Intersector1Register();
    BVH4Intersector4SingleRegister();
    BVH4Intersector4ChunkRegister();
    BVH4Intersector4HybridRegister();
#endif

#if defined(__AVX__)
    BVH4Intersector8SingleRegister();
    BVH4Intersector8ChunkRegister();
    BVH4Intersector8HybridRegister();
    BVH4Intersector1AVXRegister();
#endif
  }

  void BVH4::init(size_t numNodes, size_t numTriangles)
  {
    root = 0;
    alloc.clear();
  }

  void BVH4::print() { 
    BVHStatisticsT<BVH4>(this).print(); 
  }
}


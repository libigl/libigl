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

#include "bvh2.h"
#include "../geometry/triangles.h"
#include "../builders/heuristics.h"
#include "../builders/bvh_builder.h"
#include "../builders/bvh_statistics.h"

namespace embree
{
  /*! intersector registration functions */
  void BVH2Intersector1Register();
  void BVH2Intersector4ChunkRegister();
  void BVH2Intersector8ChunkRegister();
  
  void BVH2Register () 
  {
    /* register acceleration structure */
    TriangleMesh::accels.add("bvh2",BVH2::create);
    VirtualScene   ::accels.add("bvh2",BVH2::create);
    TriangleMesh::accels.setDefaultTriangle("bvh2","triangle4");
    VirtualScene   ::accels.setDefaultTriangle("bvh2","virtual");
    TriangleMesh::builders.setDefaultBuilder ("bvh2","objectsplit");
    VirtualScene   ::builders.setDefaultBuilder ("bvh2","objectsplit");

    /* register triangle mesh builders */
    TriangleMesh::builders.add("bvh2","triangle1i","objectsplit", BVHBuilderT<BVH2, HeuristicBinning<0> >::create,1,inf);
    TriangleMesh::builders.add("bvh2","triangle1v","objectsplit", BVHBuilderT<BVH2, HeuristicBinning<0> >::create,1,inf);
    TriangleMesh::builders.add("bvh2","triangle1", "objectsplit", BVHBuilderT<BVH2, HeuristicBinning<0> >::create,1,inf);

#if defined (__SSE__)
    TriangleMesh::builders.add("bvh2","triangle4i","objectsplit", BVHBuilderT<BVH2, HeuristicBinning<2> >::create,1,inf);
    TriangleMesh::builders.add("bvh2","triangle4v","objectsplit", BVHBuilderT<BVH2, HeuristicBinning<2> >::create,1,inf);
    TriangleMesh::builders.add("bvh2","triangle4", "objectsplit", BVHBuilderT<BVH2, HeuristicBinning<2> >::create,1,inf);
#endif

#if defined (__AVX__)
    TriangleMesh::builders.add("bvh2","triangle8", "objectsplit", BVHBuilderT<BVH2, HeuristicBinning<3> >::create,1,inf);
#endif

    TriangleMesh::builders.add("bvh2","triangle1i","spatialsplit", BVHBuilderT<BVH2, HeuristicSpatial<0> >::create,1,inf);
    TriangleMesh::builders.add("bvh2","triangle1v","spatialsplit", BVHBuilderT<BVH2, HeuristicSpatial<0> >::create,1,inf);
    TriangleMesh::builders.add("bvh2","triangle1", "spatialsplit", BVHBuilderT<BVH2, HeuristicSpatial<0> >::create,1,inf);

#if defined (__SSE__)
    TriangleMesh::builders.add("bvh2","triangle4i","spatialsplit", BVHBuilderT<BVH2, HeuristicSpatial<2> >::create,1,inf);
    TriangleMesh::builders.add("bvh2","triangle4v","spatialsplit", BVHBuilderT<BVH2, HeuristicSpatial<2> >::create,1,inf);
    TriangleMesh::builders.add("bvh2","triangle4", "spatialsplit", BVHBuilderT<BVH2, HeuristicSpatial<2> >::create,1,inf);
#endif

#if defined (__AVX__)
    TriangleMesh::builders.add("bvh2","triangle8", "spatialsplit", BVHBuilderT<BVH2, HeuristicSpatial<3> >::create,1,inf);
#endif

    /* register virtual object builders */
    VirtualScene::builders.add("bvh2","virtual", "objectsplit", BVHBuilderT<BVH2, HeuristicBinning<0> >::create,1,inf);

#if defined(__SSE__)
    BVH2Intersector1Register();
    BVH2Intersector4ChunkRegister();
#endif

#if defined (__AVX__)
    BVH2Intersector8ChunkRegister();
#endif
  }

  void BVH2::init(size_t numNodes, size_t numTriangles)
  {
    root = 0;
    alloc.clear();
  }

  void BVH2::print() { 
    BVHStatisticsT<BVH2>(this).print(); 
  }
}

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

#include "bvh4i.h"
#include "../geometry/triangles.h"
#include "../builders/heuristics.h"
#include "../builders/bvh_builder.h"
#include "../builders/bvh_statistics.h"

namespace embree
{
  /*! intersector registration functions */
  void BVH4iIntersector1Register();
  
  void BVH4iRegister () 
  {
    /* register acceleration structure */
    TriangleMesh::accels.add("bvh4i",BVH4i::create);
    VirtualScene   ::accels.add("bvh4i",BVH4i::create);
    TriangleMesh::accels.setDefaultTriangle("bvh4i","triangle4");
    VirtualScene   ::accels.setDefaultTriangle("bvh4i","virtual");
    TriangleMesh::builders.setDefaultBuilder ("bvh4i","objectsplit");
    VirtualScene   ::builders.setDefaultBuilder ("bvh4i","objectsplit");

    /* register triangle mesh builders */
    TriangleMesh::builders.add("bvh4i","triangle1i","objectsplit" ,BVHBuilderT<BVH4i, HeuristicBinning<0> >::create,1,inf);
    TriangleMesh::builders.add("bvh4i","triangle1v","objectsplit" ,BVHBuilderT<BVH4i, HeuristicBinning<0> >::create,1,inf);
    TriangleMesh::builders.add("bvh4i","triangle4i","objectsplit" ,BVHBuilderT<BVH4i, HeuristicBinning<2> >::create,1,inf);
    TriangleMesh::builders.add("bvh4i","triangle4v","objectsplit" ,BVHBuilderT<BVH4i, HeuristicBinning<2> >::create,1,inf);
    TriangleMesh::builders.add("bvh4i","triangle1" ,"objectsplit", BVHBuilderT<BVH4i, HeuristicBinning<0> >::create,1,inf);
    TriangleMesh::builders.add("bvh4i","triangle4" ,"objectsplit", BVHBuilderT<BVH4i, HeuristicBinning<2> >::create,1,inf);
#if defined (__AVX__)
    TriangleMesh::builders.add("bvh4i","triangle8" ,"objectsplit", BVHBuilderT<BVH4i, HeuristicBinning<3> >::create,1,inf);
#endif

    TriangleMesh::builders.add("bvh4i","triangle1i","spatialsplit" ,BVHBuilderT<BVH4i, HeuristicSpatial<0> >::create,1,inf);
    TriangleMesh::builders.add("bvh4i","triangle1v","spatialsplit" ,BVHBuilderT<BVH4i, HeuristicSpatial<0> >::create,1,inf);
    TriangleMesh::builders.add("bvh4i","triangle4i","spatialsplit" ,BVHBuilderT<BVH4i, HeuristicSpatial<2> >::create,1,inf);
    TriangleMesh::builders.add("bvh4i","triangle4v","spatialsplit" ,BVHBuilderT<BVH4i, HeuristicSpatial<2> >::create,1,inf);
    TriangleMesh::builders.add("bvh4i","triangle1" ,"spatialsplit", BVHBuilderT<BVH4i, HeuristicSpatial<0> >::create,1,inf);
    TriangleMesh::builders.add("bvh4i","triangle4" ,"spatialsplit", BVHBuilderT<BVH4i, HeuristicSpatial<2> >::create,1,inf);
#if defined (__AVX__)
    TriangleMesh::builders.add("bvh4i","triangle8" ,"spatialsplit", BVHBuilderT<BVH4i, HeuristicSpatial<3> >::create,1,inf);
#endif
    
    /* register virtual object builders */
    VirtualScene::builders.add("bvh4i","virtual"   ,"objectsplit", BVHBuilderT<BVH4i, HeuristicBinning<0> >::create,1,inf);

    /* register intersectors */
    BVH4iIntersector1Register();
  }

  void BVH4i::init(size_t numNodes, size_t numTriangles)
  {
    root = 0;
    alloc_nodes->init(numNodes*sizeof(Node));
    alloc_tris ->init(numTriangles*trity.bytes);
  }

  void BVH4i::print() { 
    BVHStatisticsT<BVH4i>(this).print(); 
  }
}

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

#include "sys/platform.h"

#if defined(__AVX__)

#include "bvh8.h"
#include "../geometry/triangles.h"
#include "../builders/heuristics.h"
#include "../builders/bvh_builder.h"
#include "../builders/bvh_statistics.h"

namespace embree
{
  /*! intersector registration functions */
  void BVH8Intersector1Register();
    
  void BVH8Register () 
  {
    /* register acceleration structure */
    TriangleMesh::accels.add("bvh8",BVH8::create);
    VirtualScene   ::accels.add("bvh8",BVH8::create);
    TriangleMesh::accels.setDefaultTriangle("bvh8","triangle4");
    VirtualScene   ::accels.setDefaultTriangle("bvh8","virtual");
    TriangleMesh::builders.setDefaultBuilder ("bvh8","objectsplit");
    VirtualScene   ::builders.setDefaultBuilder ("bvh8","objectsplit");

    /* register triangle mesh builders */
    TriangleMesh::builders.add("bvh8","triangle1i","objectsplit" ,BVHBuilderT<BVH8, HeuristicBinning<0> >::create,1,inf);
    TriangleMesh::builders.add("bvh8","triangle1v","objectsplit" ,BVHBuilderT<BVH8, HeuristicBinning<0> >::create,1,inf);
    TriangleMesh::builders.add("bvh8","triangle4i","objectsplit" ,BVHBuilderT<BVH8, HeuristicBinning<2> >::create,1,inf);
    TriangleMesh::builders.add("bvh8","triangle4v","objectsplit" ,BVHBuilderT<BVH8, HeuristicBinning<2> >::create,1,inf);
    TriangleMesh::builders.add("bvh8","triangle1" ,"objectsplit", BVHBuilderT<BVH8, HeuristicBinning<0> >::create,1,inf);
    TriangleMesh::builders.add("bvh8","triangle4" ,"objectsplit", BVHBuilderT<BVH8, HeuristicBinning<2> >::create,1,inf);
    TriangleMesh::builders.add("bvh8","triangle8" ,"objectsplit", BVHBuilderT<BVH8, HeuristicBinning<3> >::create,1,inf);

    TriangleMesh::builders.add("bvh8","triangle1i","spatialsplit" ,BVHBuilderT<BVH8, HeuristicSpatial<0> >::create,1,inf);
    TriangleMesh::builders.add("bvh8","triangle1v","spatialsplit" ,BVHBuilderT<BVH8, HeuristicSpatial<0> >::create,1,inf);
    TriangleMesh::builders.add("bvh8","triangle4i","spatialsplit" ,BVHBuilderT<BVH8, HeuristicSpatial<2> >::create,1,inf);
    TriangleMesh::builders.add("bvh8","triangle4v","spatialsplit" ,BVHBuilderT<BVH8, HeuristicSpatial<2> >::create,1,inf);
    TriangleMesh::builders.add("bvh8","triangle1" ,"spatialsplit", BVHBuilderT<BVH8, HeuristicSpatial<0> >::create,1,inf);
    TriangleMesh::builders.add("bvh8","triangle4" ,"spatialsplit", BVHBuilderT<BVH8, HeuristicSpatial<2> >::create,1,inf);
    TriangleMesh::builders.add("bvh8","triangle8" ,"spatialsplit", BVHBuilderT<BVH8, HeuristicSpatial<3> >::create,1,inf);

    /* register virtual object builders */
    VirtualScene::builders.add("bvh8","virtual","objectsplit", BVHBuilderT<BVH8, HeuristicBinning<0> >::create,1,inf);

    /* register intersectors */
    BVH8Intersector1Register();
  }
  
  void BVH8::init(size_t numNodes, size_t numTriangles) {
    root = NULL;
    alloc.clear();
  }

  void BVH8::print() { 
    BVHStatisticsT<BVH8>(this).print(); 
  }
}

#endif


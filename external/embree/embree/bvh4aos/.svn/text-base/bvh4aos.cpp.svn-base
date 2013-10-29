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

#include "bvh4aos.h"
#include "../geometry/triangles.h"
#include "../builders/heuristics.h"
#include "../builders/bvh_builder.h"
#include "../builders/bvh_statistics.h"

#include "builder/bvh4aos_builder.h"

namespace embree
{
  const mic_f BVH4AOS::initial_node = select(0x8888,cast_to_mic_f(mic_i(leaf_mask)),mic_f(1E10));

  /*! intersector registration functions */
  void BVH4AOSIntersector1Register();
  void BVH4AOSIntersector16ChunkRegister();
  
  void BVH4AOSTriangle1Intersector1Register();
  void BVH4AOSTriangle1Intersector16SingleRegister();
  void BVH4AOSTriangle1Intersector16ChunkRegister();
  void BVH4AOSTriangle1Intersector16HybridRegister();

  void BVH4AOSTriangle1Intersector1RefRegister();
  void BVH4AOSTriangle1Intersector16ChunkRefRegister();
  void BVH4AOSTriangle1Intersector16HybridRefRegister();
  
  void BVH4AOSRegister () 
  {
    /* register acceleration structure */
    TriangleMesh::accels.add("bvh4aos",BVH4AOS::create);
    VirtualScene   ::accels.add("bvh4aos",BVH4AOS::create);
    TriangleMesh::accels.setDefaultAccel("bvh4aos");
    VirtualScene   ::accels.setDefaultAccel("bvh4aos");
    TriangleMesh::accels.setDefaultTriangle("bvh4aos","triangle1");
    VirtualScene   ::accels.setDefaultTriangle("bvh4aos","virtual");
    //TriangleMesh::builders.setDefaultBuilder ("bvh4aos","objectsplit");
    TriangleMesh::builders.setDefaultBuilder ("bvh4aos","bowpass");
    VirtualScene   ::builders.setDefaultBuilder ("bvh4aos","objectsplit");

    /* register triangle mesh builders */
    //TriangleMesh::builders.add("bvh4aos","triangle1i","objectsplit" ,BVHBuilderT<BVH4AOS, HeuristicBinning<0> >::create,1,inf);
    //TriangleMesh::builders.add("bvh4aos","triangle1v","objectsplit" ,BVHBuilderT<BVH4AOS, HeuristicBinning<0> >::create,4,4);
    //TriangleMesh::builders.add("bvh4aos","triangle4i","objectsplit" ,BVHBuilderT<BVH4AOS, HeuristicBinning<2> >::create,1,inf);
    //TriangleMesh::builders.add("bvh4aos","triangle4v","objectsplit" ,BVHBuilderT<BVH4AOS, HeuristicBinning<2> >::create,1,inf);
    TriangleMesh::builders.add("bvh4aos","triangle1","objectsplit" ,BVHBuilderT<BVH4AOS, HeuristicBinning<0> >::create,4,4);
    //TriangleMesh::builders.add("bvh4aos","triangle4" ,"objectsplit", BVHBuilderT<BVH4AOS, HeuristicBinning<2> >::create,1,inf);
    //TriangleMesh::builders.add("bvh4aos","triangle8" ,"objectsplit", BVHBuilderT<BVH4AOS, HeuristicBinning<3> >::create,1,inf);

    //TriangleMesh::builders.add("bvh4aos","triangle1i","spatialsplit" ,BVHBuilderT<BVH4AOS, HeuristicSpatial<0> >::create,1,inf);
    //TriangleMesh::builders.add("bvh4aos","triangle1v","spatialsplit" ,BVHBuilderT<BVH4AOS, HeuristicSpatial<0> >::create,4,4);
    //TriangleMesh::builders.add("bvh4aos","triangle4i","spatialsplit" ,BVHBuilderT<BVH4AOS, HeuristicSpatial<2> >::create,1,inf);
    //TriangleMesh::builders.add("bvh4aos","triangle4v","spatialsplit" ,BVHBuilderT<BVH4AOS, HeuristicSpatial<2> >::create,1,inf);
    TriangleMesh::builders.add("bvh4aos","triangle1","spatialsplit" ,BVHBuilderT<BVH4AOS, HeuristicSpatial<0> >::create,4,4);
    //TriangleMesh::builders.add("bvh4aos","triangle4" ,"spatialsplit", BVHBuilderT<BVH4AOS, HeuristicSpatial<2> >::create,1,inf);
    //TriangleMesh::builders.add("bvh4aos","triangle8" ,"spatialsplit", BVHBuilderT<BVH4AOS, HeuristicSpatial<3> >::create,1,inf);

    TriangleMesh::builders.add("bvh4aos","triangle1","bowpass" ,BVH4AOSBuilder::create,4,4);
    TriangleMesh::builders.add("bvh4aos","triangle1","bowpass_fast" ,BVH4AOSBuilder::create_fast,4,4);

    /* register virtual object builders */
    VirtualScene::builders.add("bvh4aos","virtual","objectsplit", BVHBuilderT<BVH4AOS, HeuristicBinning<0> >::create,1,inf);

    /* register intersectors */
    BVH4AOSIntersector1Register();
    BVH4AOSIntersector16ChunkRegister();

    BVH4AOSTriangle1Intersector1Register();
    BVH4AOSTriangle1Intersector16SingleRegister();
    BVH4AOSTriangle1Intersector16ChunkRegister();
    BVH4AOSTriangle1Intersector16HybridRegister();
    
    BVH4AOSTriangle1Intersector1RefRegister();
    BVH4AOSTriangle1Intersector16ChunkRefRegister();
    BVH4AOSTriangle1Intersector16HybridRefRegister();
  }

  void BVH4AOS::init(size_t numNodes, size_t numTriangles)
  {
    root = 0;
    size_t bytesNode = max(size_t(1 << index_shift),sizeof(Node));
    size_t bytesTri  = max(size_t(1 << index_shift),trity.bytes);
    alloc_nodes->init(numNodes*bytesNode);
    alloc_tris ->init(numTriangles*bytesTri);
  }

  void BVH4AOS::resize (size_t bytesNodeArray, size_t bytesTriangleArray)
  {
    root = NULL;
    alloc_nodes->init(bytesNodeArray);
    alloc_tris ->init(bytesTriangleArray);
  }

  void BVH4AOS::print() { 
    BVHStatisticsT<BVH4AOS>(this).print(); 
  }
}

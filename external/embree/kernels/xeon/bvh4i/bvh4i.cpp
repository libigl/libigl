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

#include "geometry/triangle1.h"
#include "geometry/triangle4.h"
#include "geometry/triangle1v.h"
#include "geometry/triangle4v.h"
#include "geometry/triangle8.h"

#include "common/accelinstance.h"

namespace embree
{
 __aligned(64) Vec3fa initQBVHNode[2] = { 
    Vec3fa(pos_inf,pos_inf,pos_inf,(int)(1 << 31)),
    Vec3fa(neg_inf,neg_inf,neg_inf,(int)(1 << 31))
  };

  void BVH4iBuilderRegister ();
  void BVH4iBuilderMortonRegister ();
  void BVH4iBuilderMortonEnhancedRegister ();

  DECLARE_SYMBOL(Accel::Intersector1,BVH4iTriangle1Intersector1Moeller);
  DECLARE_SYMBOL(Accel::Intersector1,BVH4iTriangle4Intersector1Moeller);
  DECLARE_SYMBOL(Accel::Intersector1,BVH4iTriangle1vIntersector1Pluecker);
  DECLARE_SYMBOL(Accel::Intersector1,BVH4iTriangle4vIntersector1Pluecker);
  DECLARE_SYMBOL(Accel::Intersector1,BVH4iVirtualIntersector1);

  // -----
  DECLARE_SYMBOL(Accel::Intersector1,BVH4iTriangle1Intersector1ScalarMoeller);
  DECLARE_SYMBOL(Accel::Intersector1,BVH4iVirtualIntersector1Scalar);
  // -----

  DECLARE_SYMBOL(Accel::Intersector4,BVH4iTriangle1Intersector4ChunkMoeller);
  DECLARE_SYMBOL(Accel::Intersector4,BVH4iTriangle4Intersector4ChunkMoeller);
  DECLARE_SYMBOL(Accel::Intersector4,BVH4iTriangle1vIntersector4ChunkPluecker);
  DECLARE_SYMBOL(Accel::Intersector4,BVH4iTriangle4vIntersector4ChunkPluecker);
  DECLARE_SYMBOL(Accel::Intersector4,BVH4iVirtualIntersector4Chunk);

  DECLARE_SYMBOL(Accel::Intersector8,BVH4iTriangle1Intersector8ChunkMoeller);
  DECLARE_SYMBOL(Accel::Intersector8,BVH4iTriangle4Intersector8ChunkMoeller);
  DECLARE_SYMBOL(Accel::Intersector8,BVH4iTriangle1vIntersector8ChunkPluecker);
  DECLARE_SYMBOL(Accel::Intersector8,BVH4iTriangle4vIntersector8ChunkPluecker);
  DECLARE_SYMBOL(Accel::Intersector8,BVH4iVirtualIntersector8Chunk);

  DECLARE_SYMBOL(Accel::Intersector8,BVH4iTriangle4Intersector8HybridMoeller);
  DECLARE_SYMBOL(Accel::Intersector8,BVH4iTriangle8Intersector8HybridMoeller);
  DECLARE_SYMBOL(Accel::Intersector1,BVH4iTriangle8Intersector1Moeller);
  DECLARE_SYMBOL(Accel::Intersector4,BVH4iTriangle8Intersector4ChunkMoeller);
  DECLARE_SYMBOL(Accel::Intersector8,BVH4iTriangle8Intersector8ChunkMoeller);


#if defined(__TARGET_AVX2__)
  extern Accel::Intersector8 BVH4iTriangle1Intersector8ChunkAVX2;
#endif

  DECLARE_BUILDER(BVH4iTriangle1BuilderObjectSplit4Fast);
  DECLARE_BUILDER(BVH4iTriangle1BuilderMorton);
  DECLARE_BUILDER(BVH4iTriangle1BuilderMortonEnhanced);
  
  Builder* BVH4iBuilderObjectSplit1 (void* accel, BuildSource* source, void* geometry, const size_t minLeafSize, const size_t maxLeafSize);
  Builder* BVH4iBuilderObjectSplit4 (void* accel, BuildSource* source, void* geometry, const size_t minLeafSize, const size_t maxLeafSize);
  Builder* BVH4iBuilderSpatialSplit1 (void* accel, BuildSource* source, void* geometry, const size_t minLeafSize, const size_t maxLeafSize);
  Builder* BVH4iBuilderSpatialSplit4 (void* accel, BuildSource* source, void* geometry, const size_t minLeafSize, const size_t maxLeafSize);

  Builder* BVH4iBuilderObjectSplit8 (void* accel, BuildSource* source, void* geometry, const size_t minLeafSize, const size_t maxLeafSize);

  Builder* BVH4iBuilderSpatialSplit8 (void* accel, BuildSource* source, void* geometry, const size_t minLeafSize, const size_t maxLeafSize);


  void BVH4iRegister () 
  {
    int features = getCPUFeatures();

    SELECT_SYMBOL_AVX_AVX2(features,BVH4iTriangle1BuilderObjectSplit4Fast);
    SELECT_SYMBOL_AVX     (features,BVH4iTriangle1BuilderMorton);
    SELECT_SYMBOL_AVX     (features,BVH4iTriangle1BuilderMortonEnhanced);

    SELECT_SYMBOL_DEFAULT_AVX_AVX2(features,BVH4iTriangle1Intersector1Moeller);
    SELECT_SYMBOL_DEFAULT_AVX_AVX2(features,BVH4iTriangle4Intersector1Moeller);
    SELECT_SYMBOL_DEFAULT_AVX_AVX2(features,BVH4iTriangle1vIntersector1Pluecker);
    SELECT_SYMBOL_DEFAULT_AVX_AVX2(features,BVH4iTriangle4vIntersector1Pluecker);
    SELECT_SYMBOL_DEFAULT_AVX_AVX2(features,BVH4iVirtualIntersector1);

    // -----
    SELECT_SYMBOL_DEFAULT_AVX_AVX2(features,BVH4iTriangle1Intersector1ScalarMoeller);
    SELECT_SYMBOL_DEFAULT_AVX_AVX2(features,BVH4iVirtualIntersector1Scalar);
    // -----

    SELECT_SYMBOL_DEFAULT_AVX_AVX2(features,BVH4iTriangle1Intersector4ChunkMoeller);
    SELECT_SYMBOL_DEFAULT_AVX_AVX2(features,BVH4iTriangle4Intersector4ChunkMoeller);
    SELECT_SYMBOL_DEFAULT_AVX_AVX2(features,BVH4iTriangle1vIntersector4ChunkPluecker);
    SELECT_SYMBOL_DEFAULT_AVX_AVX2(features,BVH4iTriangle4vIntersector4ChunkPluecker);
    SELECT_SYMBOL_DEFAULT_AVX_AVX2(features,BVH4iVirtualIntersector4Chunk);

    SELECT_SYMBOL_AVX_AVX2(features,BVH4iTriangle1Intersector8ChunkMoeller);
    SELECT_SYMBOL_AVX_AVX2(features,BVH4iTriangle4Intersector8ChunkMoeller);
    SELECT_SYMBOL_AVX_AVX2(features,BVH4iTriangle1vIntersector8ChunkPluecker);
    SELECT_SYMBOL_AVX_AVX2(features,BVH4iTriangle4vIntersector8ChunkPluecker);
    SELECT_SYMBOL_AVX_AVX2(features,BVH4iVirtualIntersector8Chunk);

    SELECT_SYMBOL_AVX_AVX2(features,BVH4iTriangle4Intersector8HybridMoeller);
    SELECT_SYMBOL_AVX_AVX2(features,BVH4iTriangle8Intersector8HybridMoeller);
    SELECT_SYMBOL_AVX_AVX2(features,BVH4iTriangle8Intersector1Moeller);
    SELECT_SYMBOL_AVX_AVX2(features,BVH4iTriangle8Intersector4ChunkMoeller);
    SELECT_SYMBOL_AVX_AVX2(features,BVH4iTriangle8Intersector8ChunkMoeller);

  }


  Accel::Intersectors BVH4iTriangle1Intersectors(BVH4i* bvh)
  {
    Accel::Intersectors intersectors;
    intersectors.ptr = bvh;
    if (g_traverser == "scalar" ) 
      intersectors.intersector1 = BVH4iTriangle1Intersector1ScalarMoeller;
    else
      intersectors.intersector1 = BVH4iTriangle1Intersector1Moeller;
    intersectors.intersector4 = BVH4iTriangle1Intersector4ChunkMoeller;
    intersectors.intersector8 = BVH4iTriangle1Intersector8ChunkMoeller;
    intersectors.intersector16 = NULL;
    return intersectors;
  }

  Accel::Intersectors BVH4iTriangle4Intersectors(BVH4i* bvh)
  {
    Accel::Intersectors intersectors;
    intersectors.ptr = bvh;
    intersectors.intersector1 = BVH4iTriangle4Intersector1Moeller;
    intersectors.intersector4 = BVH4iTriangle4Intersector4ChunkMoeller;
    intersectors.intersector8 = BVH4iTriangle4Intersector8HybridMoeller;
    intersectors.intersector16 = NULL;
    return intersectors;
  }

  Accel::Intersectors BVH4iTriangle1vIntersectors(BVH4i* bvh)
  {
    Accel::Intersectors intersectors;
    intersectors.ptr = bvh;
    intersectors.intersector1 = BVH4iTriangle1vIntersector1Pluecker;
    intersectors.intersector4 = BVH4iTriangle1vIntersector4ChunkPluecker;
    intersectors.intersector8 = BVH4iTriangle1vIntersector8ChunkPluecker;
    intersectors.intersector16 = NULL;
    return intersectors;
  }

  Accel::Intersectors BVH4iTriangle4vIntersectors(BVH4i* bvh)
  {
    Accel::Intersectors intersectors;
    intersectors.ptr = bvh;
    intersectors.intersector1 = BVH4iTriangle4vIntersector1Pluecker;
    intersectors.intersector4 = BVH4iTriangle4vIntersector4ChunkPluecker;
    intersectors.intersector8 = BVH4iTriangle4vIntersector8ChunkPluecker;
    intersectors.intersector16 = NULL;
    return intersectors;
  }

  Accel::Intersectors BVH4iTriangle8IntersectorsHybrid(BVH4i* bvh)
  {
    Accel::Intersectors intersectors;
    intersectors.ptr = bvh;
    intersectors.intersector1 = BVH4iTriangle8Intersector1Moeller;
    intersectors.intersector4 = BVH4iTriangle8Intersector4ChunkMoeller;
    intersectors.intersector8 = BVH4iTriangle8Intersector8HybridMoeller;
    intersectors.intersector16 = NULL;
    return intersectors;
  }

#if defined (__TARGET_AVX__)

  Accel* BVH4i::BVH4iTriangle8(Scene* scene)
  { 
    BVH4i* accel = new BVH4i(SceneTriangle8::type,scene);

    Accel::Intersectors intersectors;
    if      (g_traverser == "default") intersectors = BVH4iTriangle8IntersectorsHybrid(accel);
    else if (g_traverser == "hybrid" ) intersectors = BVH4iTriangle8IntersectorsHybrid(accel);
    else throw std::runtime_error("unknown traverser "+g_traverser+" for BVH4i<Triangle8>");
   
    Builder* builder = NULL;
    if      (g_builder == "default"     ) builder = BVH4iBuilderObjectSplit8(accel,&scene->flat_triangle_source_1,scene,1,inf);
    else if (g_builder == "spatialsplit") builder = BVH4iBuilderSpatialSplit8(accel,&scene->flat_triangle_source_1,scene,1,inf);
    else if (g_builder == "objectsplit" ) builder = BVH4iBuilderObjectSplit8(accel,&scene->flat_triangle_source_1,scene,1,inf);
    else throw std::runtime_error("unknown builder "+g_builder+" for BVH4i<Triangle8>");

    return new AccelInstance(accel,builder,intersectors);
  }
#endif


  Accel* BVH4i::BVH4iTriangle1(Scene* scene)
  { 
    BVH4i* accel = new BVH4i(SceneTriangle1::type);
    
    Builder* builder = NULL;
    if      (g_builder == "default"     ) builder = BVH4iBuilderObjectSplit1(accel,&scene->flat_triangle_source_1,scene,1,inf);
    else if (g_builder == "spatialsplit") builder = BVH4iBuilderSpatialSplit1(accel,&scene->flat_triangle_source_1,scene,1,inf);
    else if (g_builder == "objectsplit" ) builder = BVH4iBuilderObjectSplit1(accel,&scene->flat_triangle_source_1,scene,1,inf);
    else if (g_builder == "objectsplit_fast") builder = BVH4iTriangle1BuilderObjectSplit4Fast(accel,&scene->flat_triangle_source_1,scene,1,inf);
    else if (g_builder == "morton"          ) builder = BVH4iTriangle1BuilderMorton(accel,&scene->flat_triangle_source_1,scene,1,inf);
    else if (g_builder == "morton.enhanced" ) builder = BVH4iTriangle1BuilderMortonEnhanced(accel,&scene->flat_triangle_source_1,scene,1,inf);
    else throw std::runtime_error("unknown builder "+g_builder+" for BVH4i<Triangle1>");
    
    Accel::Intersectors intersectors = BVH4iTriangle1Intersectors(accel);
    return new AccelInstance(accel,builder,intersectors);
  }
  
  Accel* BVH4i::BVH4iTriangle4(Scene* scene)
  { 
    BVH4i* accel = new BVH4i(SceneTriangle4::type);
    
    Builder* builder = NULL;
    if      (g_builder == "default"     ) builder = BVH4iBuilderObjectSplit4(accel,&scene->flat_triangle_source_1,scene,1,inf);
    else if (g_builder == "spatialsplit") builder = BVH4iBuilderSpatialSplit4(accel,&scene->flat_triangle_source_1,scene,1,inf);
    else if (g_builder == "objectsplit" ) builder = BVH4iBuilderObjectSplit4(accel,&scene->flat_triangle_source_1,scene,1,inf);
    else throw std::runtime_error("unknown builder "+g_builder+" for BVH4i<Triangle4>");
    
    Accel::Intersectors intersectors = BVH4iTriangle4Intersectors(accel);
    return new AccelInstance(accel,builder,intersectors);
  }
  
  Accel* BVH4i::BVH4iTriangle1_v1(Scene* scene)
  { 
    BVH4i* accel = new BVH4i(SceneTriangle1::type);
    Builder* builder = BVH4iTriangle1BuilderObjectSplit4Fast(accel,&scene->flat_triangle_source_1,scene,1,inf);
    
    Accel::Intersectors intersectors = BVH4iTriangle1Intersectors(accel);
    return new AccelInstance(accel,builder,intersectors);
  }
  
  Accel* BVH4i::BVH4iTriangle1_v2(Scene* scene)
  { 
    BVH4i* accel = new BVH4i(SceneTriangle1::type);
    Builder* builder = BVH4iTriangle1BuilderObjectSplit4Fast(accel,&scene->flat_triangle_source_1,scene,1,inf);
    
    Accel::Intersectors intersectors = BVH4iTriangle1Intersectors(accel);
#if defined(__TARGET_AVX2__)
    intersectors.intersector8 = BVH4iTriangle1Intersector8ChunkAVX2;
#endif
    return new AccelInstance(accel,builder,intersectors);
  }
  
  Accel* BVH4i::BVH4iTriangle1_morton(Scene* scene)
  { 
    BVH4i* accel = new BVH4i(SceneTriangle1::type);
    Builder* builder = BVH4iTriangle1BuilderMorton(accel,&scene->flat_triangle_source_1,scene,1,inf);
    Accel::Intersectors intersectors = BVH4iTriangle1Intersectors(accel);
    return new AccelInstance(accel,builder,intersectors);
  }
  
  Accel* BVH4i::BVH4iTriangle1_morton_enhanced(Scene* scene)
  { 
    BVH4i* accel = new BVH4i(SceneTriangle1::type);
    Builder* builder = BVH4iTriangle1BuilderMortonEnhanced(accel,&scene->flat_triangle_source_1,scene,1,inf);
    Accel::Intersectors intersectors = BVH4iTriangle1Intersectors(accel);
    return new AccelInstance(accel,builder,intersectors);
  }
  
  Accel* BVH4i::BVH4iTriangle1(TriangleMeshScene::TriangleMesh* mesh)
  {
    BVH4i* accel = new BVH4i(TriangleMeshTriangle1::type);

    Builder* builder = NULL;
    if      (g_builder == "default"     ) builder = BVH4iBuilderObjectSplit1(accel,mesh,mesh,1,inf);
    else if (g_builder == "spatialsplit") builder = BVH4iBuilderSpatialSplit1(accel,mesh,mesh,1,inf);
    else if (g_builder == "objectsplit" ) builder = BVH4iBuilderObjectSplit1(accel,mesh,mesh,1,inf);
    else throw std::runtime_error("unknown builder "+g_builder+" for BVH4i<Triangle1>");

    Accel::Intersectors intersectors = BVH4iTriangle1Intersectors(accel);
    return new AccelInstance(accel,builder,intersectors);
  }

  Accel* BVH4i::BVH4iTriangle4(TriangleMeshScene::TriangleMesh* mesh)
  {
    BVH4i* accel = new BVH4i(TriangleMeshTriangle4::type);

    Builder* builder = NULL;
    if      (g_builder == "default"     ) builder = BVH4iBuilderObjectSplit4(accel,mesh,mesh,1,inf);
    else if (g_builder == "spatialsplit") builder = BVH4iBuilderSpatialSplit4(accel,mesh,mesh,1,inf);
    else if (g_builder == "objectsplit" ) builder = BVH4iBuilderObjectSplit4(accel,mesh,mesh,1,inf);
    else throw std::runtime_error("unknown builder "+g_builder+" for BVH4i<Triangle4>");

    Accel::Intersectors intersectors = BVH4iTriangle4Intersectors(accel);
    return new AccelInstance(accel,builder,intersectors);
  }

  Accel* BVH4i::BVH4iTriangle1v(TriangleMeshScene::TriangleMesh* mesh)
  {
    BVH4i* accel = new BVH4i(TriangleMeshTriangle1v::type);

    Builder* builder = NULL;
    if      (g_builder == "default"     ) builder = BVH4iBuilderObjectSplit1(accel,mesh,mesh,1,inf);
    else if (g_builder == "spatialsplit") builder = BVH4iBuilderSpatialSplit1(accel,mesh,mesh,1,inf);
    else if (g_builder == "objectsplit" ) builder = BVH4iBuilderObjectSplit1(accel,mesh,mesh,1,inf);
    else throw std::runtime_error("unknown builder "+g_builder+" for BVH4i<Triangle1v>");

    Accel::Intersectors intersectors = BVH4iTriangle1vIntersectors(accel);
    return new AccelInstance(accel,builder,intersectors);
  }

  Accel* BVH4i::BVH4iTriangle4v(TriangleMeshScene::TriangleMesh* mesh)
  {
    BVH4i* accel = new BVH4i(TriangleMeshTriangle4v::type);

    Builder* builder = NULL;
    if      (g_builder == "default"     ) builder = BVH4iBuilderObjectSplit4(accel,mesh,mesh,1,inf);
    else if (g_builder == "spatialsplit") builder = BVH4iBuilderSpatialSplit4(accel,mesh,mesh,1,inf);
    else if (g_builder == "objectsplit" ) builder = BVH4iBuilderObjectSplit4(accel,mesh,mesh,1,inf);
    else throw std::runtime_error("unknown builder "+g_builder+" for BVH4i<Triangle4v>");
    
    Accel::Intersectors intersectors = BVH4iTriangle4vIntersectors(accel);
    return new AccelInstance(accel,builder,intersectors);
  }

  void BVH4i::init(size_t numNodes, size_t numPrimitives)
  {
    root = emptyNode;
    alloc_nodes->init(numNodes*sizeof(BVH4i::Node));
    alloc_tris ->init(numPrimitives*primTy.bytes);
  }

  void BVH4i::clearBarrier(NodeRef& node)
  {
    if (node.isBarrier()) 
      node.clearBarrier();
    else if (!node.isLeaf()) {
      Node* n = node.node(nodePtr());
      for (size_t c=0; c<4; c++)
        clearBarrier(n->child(c));
    }
  }

  float BVH4i::sah () {
    return sah(root,bounds)/area(bounds);
  }

  float BVH4i::sah (NodeRef& node, const BBox3f& bounds)
  {
    float f = bounds.empty() ? 0.0f : area(bounds);

    if (node.isNode()) 
    {
      Node* n = node.node(nodePtr());
      for (size_t c=0; c<4; c++) 
        f += sah(n->child(c),n->bounds(c));
      return f;
    }
    else 
    {
      size_t num; node.leaf(triPtr(),num);
      return f*num;
    }
  }
  
}

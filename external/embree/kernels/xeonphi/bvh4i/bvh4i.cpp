// ======================================================================== //
// Copyright 2009-2014 Intel Corporation                                    //
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
#include "bvh4i_builder.h"
#include "bvh4i_builder_morton.h"

#include "common/accelinstance.h"
#include "geometry/triangle1.h"

namespace embree
{

  __aligned(64) BVH4i::Helper BVH4i::initQBVHNode[4] = { 
    { pos_inf, pos_inf, pos_inf,BVH4i::invalidNode},
    { neg_inf, neg_inf, neg_inf,BVH4i::invalidNode},
    { pos_inf, pos_inf, pos_inf,BVH4i::invalidNode},
    { neg_inf, neg_inf, neg_inf,BVH4i::invalidNode}
   };

  /*! intersector registration functions */
  DECLARE_SYMBOL(Accel::Intersector1 ,BVH4iTriangle1Intersector1);
  DECLARE_SYMBOL(Accel::Intersector1 ,BVH4iTriangle1Intersector1NoFilter);

  DECLARE_SYMBOL(Accel::Intersector16,BVH4iTriangle1Intersector16ChunkMoeller);
  DECLARE_SYMBOL(Accel::Intersector16,BVH4iTriangle1Intersector16ChunkMoellerNoFilter);

  DECLARE_SYMBOL(Accel::Intersector16,BVH4iTriangle1Intersector16SingleMoeller);
  DECLARE_SYMBOL(Accel::Intersector16,BVH4iTriangle1Intersector16SingleMoellerNoFilter);

  DECLARE_SYMBOL(Accel::Intersector16,BVH4iTriangle1Intersector16HybridMoeller);
  DECLARE_SYMBOL(Accel::Intersector16,BVH4iTriangle1Intersector16HybridMoellerNoFilter);

  DECLARE_SYMBOL(Accel::Intersector16,BVH4iTriangle1Intersector16TestMoeller);
  DECLARE_SYMBOL(Accel::Intersector16,BVH4iTriangle1Intersector16TestMoellerNoFilter);

  DECLARE_SYMBOL(Accel::Intersector1 ,BVH4iVirtualGeometryIntersector1);
  DECLARE_SYMBOL(Accel::Intersector1 ,BVH4iVirtualGeometryIntersector1NoFilter);

  DECLARE_SYMBOL(Accel::Intersector16,BVH4iVirtualGeometryIntersector16);
  DECLARE_SYMBOL(Accel::Intersector16,BVH4iVirtualGeometryIntersector16NoFilter);


  DECLARE_SYMBOL(Accel::Intersector1 ,BVH4iTriangle1mcIntersector1);
  DECLARE_SYMBOL(Accel::Intersector1 ,BVH4iTriangle1mcIntersector1NoFilter);

  DECLARE_SYMBOL(Accel::Intersector16,BVH4iTriangle1mcIntersector16SingleMoeller);
  DECLARE_SYMBOL(Accel::Intersector16,BVH4iTriangle1mcIntersector16SingleMoellerNoFilter);

  DECLARE_SYMBOL(Accel::Intersector16,BVH4iTriangle1mcIntersector16ChunkMoeller);
  DECLARE_SYMBOL(Accel::Intersector16,BVH4iTriangle1mcIntersector16ChunkMoellerNoFilter);

  DECLARE_SYMBOL(Accel::Intersector16,BVH4iTriangle1mcIntersector16HybridMoeller);
  DECLARE_SYMBOL(Accel::Intersector16,BVH4iTriangle1mcIntersector16HybridMoellerNoFilter);

  DECLARE_SYMBOL(Accel::Intersector1 ,BVH4iSubdivMeshIntersector1);
  DECLARE_SYMBOL(Accel::Intersector1 ,BVH4iSubdivMeshIntersector1NoFilter);

  DECLARE_SYMBOL(Accel::Intersector16,BVH4iSubdivMeshIntersector16);
  DECLARE_SYMBOL(Accel::Intersector16,BVH4iSubdivMeshIntersector16NoFilter);

  void BVH4iRegister () 
  {
    int features = getCPUFeatures();

    /* default target */
    SELECT_SYMBOL_KNC(features,BVH4iTriangle1Intersector1);
    SELECT_SYMBOL_KNC(features,BVH4iTriangle1Intersector1NoFilter);

    SELECT_SYMBOL_KNC(features,BVH4iTriangle1Intersector16ChunkMoeller);
    SELECT_SYMBOL_KNC(features,BVH4iTriangle1Intersector16ChunkMoellerNoFilter);

    SELECT_SYMBOL_KNC(features,BVH4iTriangle1Intersector16SingleMoeller);
    SELECT_SYMBOL_KNC(features,BVH4iTriangle1Intersector16SingleMoellerNoFilter);

    SELECT_SYMBOL_KNC(features,BVH4iTriangle1Intersector16HybridMoeller);
    SELECT_SYMBOL_KNC(features,BVH4iTriangle1Intersector16HybridMoellerNoFilter);

    SELECT_SYMBOL_KNC(features,BVH4iTriangle1Intersector16TestMoeller);
    SELECT_SYMBOL_KNC(features,BVH4iTriangle1Intersector16TestMoellerNoFilter);

    SELECT_SYMBOL_KNC(features,BVH4iVirtualGeometryIntersector1);
    SELECT_SYMBOL_KNC(features,BVH4iVirtualGeometryIntersector1NoFilter);

    SELECT_SYMBOL_KNC(features,BVH4iVirtualGeometryIntersector16);
    SELECT_SYMBOL_KNC(features,BVH4iVirtualGeometryIntersector16NoFilter);

    SELECT_SYMBOL_KNC(features,BVH4iTriangle1mcIntersector1);
    SELECT_SYMBOL_KNC(features,BVH4iTriangle1mcIntersector1NoFilter);

    SELECT_SYMBOL_KNC(features,BVH4iTriangle1mcIntersector16SingleMoeller);
    SELECT_SYMBOL_KNC(features,BVH4iTriangle1mcIntersector16SingleMoellerNoFilter);

    SELECT_SYMBOL_KNC(features,BVH4iTriangle1mcIntersector16ChunkMoeller);
    SELECT_SYMBOL_KNC(features,BVH4iTriangle1mcIntersector16ChunkMoellerNoFilter);

    SELECT_SYMBOL_KNC(features,BVH4iTriangle1mcIntersector16HybridMoeller);
    SELECT_SYMBOL_KNC(features,BVH4iTriangle1mcIntersector16HybridMoellerNoFilter);

    SELECT_SYMBOL_KNC(features,BVH4iSubdivMeshIntersector1);
    SELECT_SYMBOL_KNC(features,BVH4iSubdivMeshIntersector1NoFilter);

    SELECT_SYMBOL_KNC(features,BVH4iSubdivMeshIntersector16);
    SELECT_SYMBOL_KNC(features,BVH4iSubdivMeshIntersector16NoFilter);

  }


  Accel::Intersectors BVH4iTriangle1Intersectors(BVH4i* bvh)
  {

    Accel::Intersectors intersectors;
    intersectors.ptr = bvh;
    intersectors.intersector1  = BVH4iTriangle1Intersector1;
    if      (g_tri_traverser == "default" || g_tri_traverser == "hybrid") 
      {
	intersectors.intersector16          = BVH4iTriangle1Intersector16HybridMoeller;
	intersectors.intersector16_filter   = BVH4iTriangle1Intersector16HybridMoeller;
	intersectors.intersector16_nofilter = BVH4iTriangle1Intersector16HybridMoellerNoFilter;
      }
    else if (g_tri_traverser == "chunk"  ) 
      {
	intersectors.intersector16          = BVH4iTriangle1Intersector16ChunkMoeller;
	intersectors.intersector16_filter   = BVH4iTriangle1Intersector16ChunkMoeller;
	intersectors.intersector16_nofilter = BVH4iTriangle1Intersector16ChunkMoellerNoFilter;

      }
    else if (g_tri_traverser == "single" ) 
      {
	intersectors.intersector16          = BVH4iTriangle1Intersector16SingleMoeller;
	intersectors.intersector16_filter   = BVH4iTriangle1Intersector16SingleMoeller;
	intersectors.intersector16_nofilter = BVH4iTriangle1Intersector16SingleMoellerNoFilter;

      }
    else if (g_tri_traverser == "test" ) 
      {
	intersectors.intersector16          = BVH4iTriangle1Intersector16TestMoeller;
	intersectors.intersector16_filter   = BVH4iTriangle1Intersector16TestMoeller;
	intersectors.intersector16_nofilter = BVH4iTriangle1Intersector16TestMoellerNoFilter;

      }
    else THROW_RUNTIME_ERROR("unknown traverser "+g_tri_traverser+" for BVH4i<Triangle1>");      
    return intersectors;
  }

  Accel::Intersectors BVH4iTriangle1mcIntersectors(BVH4i* bvh)
  {
    Accel::Intersectors intersectors;
    intersectors.ptr = bvh;
    intersectors.intersector1  = BVH4iTriangle1mcIntersector1; 
    if      (g_tri_traverser == "default" || g_tri_traverser == "hybrid") 
      {
	intersectors.intersector16          = BVH4iTriangle1mcIntersector16HybridMoeller;
	intersectors.intersector16_filter   = BVH4iTriangle1mcIntersector16HybridMoeller;
	intersectors.intersector16_nofilter = BVH4iTriangle1mcIntersector16HybridMoellerNoFilter;

      }
    else if (g_tri_traverser == "chunk"  ) 
      {
	intersectors.intersector16          = BVH4iTriangle1mcIntersector16ChunkMoeller;
	intersectors.intersector16_filter   = BVH4iTriangle1mcIntersector16ChunkMoeller;
	intersectors.intersector16_nofilter = BVH4iTriangle1mcIntersector16ChunkMoellerNoFilter;

      }
    else if (g_tri_traverser == "single" ) 
      {
	intersectors.intersector16          = BVH4iTriangle1mcIntersector16SingleMoeller;
	intersectors.intersector16_filter   = BVH4iTriangle1mcIntersector16SingleMoeller;
	intersectors.intersector16_nofilter = BVH4iTriangle1mcIntersector16SingleMoellerNoFilter;

      }
    else THROW_RUNTIME_ERROR("unknown traverser "+g_tri_traverser+" for BVH4i<Triangle1>");      
    return intersectors;
  }


  Accel::Intersectors BVH4iVirtualGeometryIntersectors(BVH4i* bvh)
  {
    Accel::Intersectors intersectors;
    intersectors.ptr = bvh;
    intersectors.intersector1  = BVH4iVirtualGeometryIntersector1;
    intersectors.intersector16 = BVH4iVirtualGeometryIntersector16;
    return intersectors;
  }

  Accel::Intersectors BVH4iSubdivMeshIntersectors(BVH4i* bvh)
  {
    Accel::Intersectors intersectors;
    intersectors.ptr = bvh;
    intersectors.intersector1  = BVH4iSubdivMeshIntersector1;
    intersectors.intersector16 = BVH4iSubdivMeshIntersector16;
    return intersectors;
  }


  Accel* BVH4i::BVH4iTriangle1ObjectSplitBinnedSAH(Scene* scene)
  { 
    BVH4i* accel = new BVH4i(SceneTriangle1::type,scene);   
    Builder* builder = BVH4iBuilder::create(accel,scene,BVH4iBuilder::BVH4I_BUILDER_DEFAULT);    
    Accel::Intersectors intersectors = BVH4iTriangle1Intersectors(accel);
    return new AccelInstance(accel,builder,intersectors);
  }

  Accel* BVH4i::BVH4iTriangle1ObjectSplitMorton(Scene* scene)
  { 
    BVH4i* accel = new BVH4i(SceneTriangle1::type,scene);   
    Builder* builder = BVH4iBuilderMorton::create(accel,scene,false);  
    Accel::Intersectors intersectors = BVH4iTriangle1Intersectors(accel);
    return new AccelInstance(accel,builder,intersectors);
  }

  Accel* BVH4i::BVH4iTriangle1ObjectSplitMorton64Bit(Scene* scene)
  { 
    BVH4i* accel = new BVH4i(SceneTriangle1::type,scene);   
    Builder* builder = BVH4iBuilderMorton64Bit::create(accel,scene);  
    Accel::Intersectors intersectors = BVH4iTriangle1Intersectors(accel);
    return new AccelInstance(accel,builder,intersectors);
  }

  Accel* BVH4i::BVH4iTriangle1ObjectSplitEnhancedMorton(Scene* scene)
  { 
    BVH4i* accel = new BVH4i(SceneTriangle1::type,scene);
    Builder* builder = BVH4iBuilderMorton::create(accel,scene,true);  
    
    Accel::Intersectors intersectors = BVH4iTriangle1Intersectors(accel);
    return new AccelInstance(accel,builder,intersectors);
  }

  Accel* BVH4i::BVH4iTriangle1PreSplitsBinnedSAH(Scene* scene)
  {
    BVH4i* accel = new BVH4i(SceneTriangle1::type,scene);
    
    Builder* builder = BVH4iBuilder::create(accel,scene,BVH4iBuilder::BVH4I_BUILDER_PRESPLITS);
    
    Accel::Intersectors intersectors = BVH4iTriangle1Intersectors(accel);
    return new AccelInstance(accel,builder,intersectors);    
  }


  Accel* BVH4i::BVH4iTriangle1MemoryConservativeBinnedSAH(Scene* scene)
  {
    BVH4i* accel = new BVH4i(SceneTriangle1::type,scene);
    
    Builder* builder = BVH4iBuilder::create(accel,scene,BVH4iBuilder::BVH4I_BUILDER_MEMORY_CONSERVATIVE);       
    Accel::Intersectors intersectors = BVH4iTriangle1mcIntersectors(accel);
    scene->needVertices = true;

    return new AccelInstance(accel,builder,intersectors);    
  }

  Accel* BVH4i::BVH4iVirtualGeometryBinnedSAH(Scene* scene)
  {
    BVH4i* accel = new BVH4i(SceneTriangle1::type,scene);    
    Builder* builder = BVH4iBuilder::create(accel,scene,BVH4iBuilder::BVH4I_BUILDER_VIRTUAL_GEOMETRY);   
    Accel::Intersectors intersectors = BVH4iVirtualGeometryIntersectors(accel);
    return new AccelInstance(accel,builder,intersectors);    
  }

  Accel* BVH4i::BVH4iSubdivMeshBinnedSAH(Scene* scene)
  {
    BVH4i* accel = new BVH4i(SceneTriangle1::type,scene);    
    Builder* builder = BVH4iBuilder::create(accel,scene,BVH4iBuilder::BVH4I_BUILDER_SUBDIV_MESH);   
    Accel::Intersectors intersectors = BVH4iSubdivMeshIntersectors(accel);
    scene->needVertices = true;
    return new AccelInstance(accel,builder,intersectors);        
  }

  BVH4i::~BVH4i()
  {
    if (qbvh)  os_free(qbvh,size_node);
    if (accel) os_free(accel,size_accel);
  }


  float BVH4i::sah () {
    return sah(root,bounds)/area(bounds);
  }

  float BVH4i::sah (NodeRef& node, BBox3fa bounds)
  {
    float f = bounds.empty() ? 0.0f : area(bounds);

    if (node.isNode()) 
    {
      Node* n = node.node(nodePtr());
      for (size_t c=0; c<BVH4i::N; c++) 
	if (n->child(c) != BVH4i::invalidNode)
	  f += sah(n->child(c),n->bounds(c));
      return f;
    }
    else 
    {
      unsigned int num; node.leaf(triPtr(),num);
      return f*num;
    }
  }

}

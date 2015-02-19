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

#include "bvh8.h"
#include "geometry/triangle4.h"
#include "geometry/triangle8.h"
#include "common/accelinstance.h"

namespace embree
{
  DECLARE_SYMBOL(Accel::Intersector1,BVH8Triangle4Intersector1Moeller);
  DECLARE_SYMBOL(Accel::Intersector4,BVH8Triangle4Intersector4HybridMoeller);
  DECLARE_SYMBOL(Accel::Intersector4,BVH8Triangle4Intersector4HybridMoellerNoFilter);
  DECLARE_SYMBOL(Accel::Intersector8,BVH8Triangle4Intersector8ChunkMoeller);
  DECLARE_SYMBOL(Accel::Intersector8,BVH8Triangle4Intersector8HybridMoeller);
  DECLARE_SYMBOL(Accel::Intersector8,BVH8Triangle4Intersector8HybridMoellerNoFilter);

  DECLARE_SYMBOL(Accel::Intersector1,BVH8Triangle8Intersector1Moeller);
  DECLARE_SYMBOL(Accel::Intersector4,BVH8Triangle8Intersector4HybridMoeller);
  DECLARE_SYMBOL(Accel::Intersector4,BVH8Triangle8Intersector4HybridMoellerNoFilter);
  DECLARE_SYMBOL(Accel::Intersector8,BVH8Triangle8Intersector8ChunkMoeller);
  DECLARE_SYMBOL(Accel::Intersector8,BVH8Triangle8Intersector8HybridMoeller);
  DECLARE_SYMBOL(Accel::Intersector8,BVH8Triangle8Intersector8HybridMoellerNoFilter);

  DECLARE_SCENE_BUILDER(BVH8Triangle4Builder);
  DECLARE_SCENE_BUILDER(BVH8Triangle8Builder);

  void BVH8Register () 
  {
    int features = getCPUFeatures();

    /* select builders */
    SELECT_SYMBOL_AVX(features,BVH8Triangle4Builder);
    SELECT_SYMBOL_AVX(features,BVH8Triangle8Builder);
 
    /* select intersectors1 */
    SELECT_SYMBOL_AVX_AVX2(features,BVH8Triangle4Intersector1Moeller);
    SELECT_SYMBOL_AVX_AVX2(features,BVH8Triangle8Intersector1Moeller);

    /* select intersectors4 */
    SELECT_SYMBOL_AVX_AVX2(features,BVH8Triangle4Intersector4HybridMoeller);
    SELECT_SYMBOL_AVX_AVX2(features,BVH8Triangle4Intersector4HybridMoellerNoFilter);
    SELECT_SYMBOL_AVX_AVX2(features,BVH8Triangle8Intersector4HybridMoeller);
    SELECT_SYMBOL_AVX_AVX2(features,BVH8Triangle8Intersector4HybridMoellerNoFilter);

    /* select intersectors8 */
    SELECT_SYMBOL_AVX_AVX2(features,BVH8Triangle4Intersector8ChunkMoeller);
    SELECT_SYMBOL_AVX_AVX2(features,BVH8Triangle4Intersector8HybridMoeller);
    SELECT_SYMBOL_AVX_AVX2(features,BVH8Triangle4Intersector8HybridMoellerNoFilter);
    SELECT_SYMBOL_AVX_AVX2(features,BVH8Triangle8Intersector8ChunkMoeller);
    SELECT_SYMBOL_AVX_AVX2(features,BVH8Triangle8Intersector8HybridMoeller);
    SELECT_SYMBOL_AVX_AVX2(features,BVH8Triangle8Intersector8HybridMoellerNoFilter);
  }

  BVH8::BVH8 (const PrimitiveType& primTy, Scene* scene)
    : primTy(primTy), scene(scene), root(emptyNode),
      numPrimitives(0), numVertices(0) {}

  BVH8::~BVH8 () {
    for (size_t i=0; i<objects.size(); i++) 
      delete objects[i];
  }

  void BVH8::init(size_t nodeSize, size_t numPrimitives, size_t numThreads)
  {
    /* allocate as much memory as likely needed and reserve conservative amounts of memory */
    size_t blockSize = LinearAllocatorPerThread::allocBlockSize;
    
    size_t numPrimBlocks = primTy.blocks(numPrimitives);
    size_t numAllocatedNodes      = min(size_t(0.6*numPrimBlocks),numPrimitives);
    size_t numAllocatedPrimitives = min(size_t(1.2*numPrimBlocks),numPrimitives);
#if defined(__X86_64__)
    size_t numReservedNodes = 2*numPrimitives;
    size_t numReservedPrimitives = 2*numPrimitives;
#else
    size_t numReservedNodes = 1.5*numAllocatedNodes;
    size_t numReservedPrimitives = 1.5*numAllocatedPrimitives;
#endif
    
    size_t bytesAllocated = numAllocatedNodes * nodeSize + numAllocatedPrimitives * primTy.bytes; // required also for parallel split stage in BVH4BuilderFast
    size_t bytesReserved  = numReservedNodes  * nodeSize + numReservedPrimitives  * primTy.bytes;
    if (numPrimitives) bytesReserved = (bytesReserved+blockSize-1)/blockSize*blockSize + numThreads*blockSize*2;

    root = emptyNode;
    bounds = empty;
    alloc.init(bytesAllocated,bytesReserved);
  }

  void BVH8::clearBarrier(NodeRef& node)
  {
    if (node.isBarrier())
      node.clearBarrier();
    else if (!node.isLeaf()) {
      Node* n = node.node();
      for (size_t c=0; c<N; c++)
        clearBarrier(n->child(c));
    }
  }

  Accel::Intersectors BVH8Triangle4Intersectors(BVH8* bvh)
  {
    Accel::Intersectors intersectors;
    intersectors.ptr = bvh;
    intersectors.intersector1 = BVH8Triangle4Intersector1Moeller;
    intersectors.intersector4_filter = BVH8Triangle4Intersector4HybridMoeller;
    intersectors.intersector4_nofilter = BVH8Triangle4Intersector4HybridMoellerNoFilter;
    intersectors.intersector8_filter = BVH8Triangle4Intersector8HybridMoeller;
    intersectors.intersector8_nofilter = BVH8Triangle4Intersector8HybridMoellerNoFilter;
    intersectors.intersector16 = NULL;
    return intersectors;
  }

  Accel::Intersectors BVH8Triangle8Intersectors(BVH8* bvh)
  {
    Accel::Intersectors intersectors;
    intersectors.ptr = bvh;
    intersectors.intersector1 = BVH8Triangle8Intersector1Moeller;
    intersectors.intersector4_filter = BVH8Triangle8Intersector4HybridMoeller;
    intersectors.intersector4_nofilter = BVH8Triangle8Intersector4HybridMoellerNoFilter;
    intersectors.intersector8_filter = BVH8Triangle8Intersector8HybridMoeller;
    intersectors.intersector8_nofilter = BVH8Triangle8Intersector8HybridMoellerNoFilter;
    intersectors.intersector16 = NULL;
    return intersectors;
  }

  Accel* BVH8::BVH8Triangle4(Scene* scene)
  { 
    BVH8* accel = new BVH8(Triangle4Type::type,scene);
    Accel::Intersectors intersectors= BVH8Triangle4Intersectors(accel);
    
    Builder* builder = NULL;
    if      (g_tri_builder == "default"     ) builder = BVH8Triangle4Builder(accel,scene,0);
    else if (g_tri_builder == "spatialsplit") builder = BVH8Triangle4Builder(accel,scene,MODE_HIGH_QUALITY);
    else if (g_tri_builder == "objectsplit" ) builder = BVH8Triangle4Builder(accel,scene,0);
    else THROW_RUNTIME_ERROR("unknown builder "+g_tri_builder+" for BVH8<Triangle4>");
    
    return new AccelInstance(accel,builder,intersectors);
  }

  Accel* BVH8::BVH8Triangle4ObjectSplit(Scene* scene)
  {
    BVH8* accel = new BVH8(Triangle4Type::type,scene);
    Accel::Intersectors intersectors= BVH8Triangle4Intersectors(accel);
    Builder* builder = BVH8Triangle4Builder(accel,scene,0);
    return new AccelInstance(accel,builder,intersectors);
  }

  Accel* BVH8::BVH8Triangle4SpatialSplit(Scene* scene)
  {
    BVH8* accel = new BVH8(Triangle4Type::type,scene);
    Accel::Intersectors intersectors= BVH8Triangle4Intersectors(accel);
    Builder* builder = BVH8Triangle4Builder(accel,scene,MODE_HIGH_QUALITY);
    return new AccelInstance(accel,builder,intersectors);
  }

  Accel* BVH8::BVH8Triangle8(Scene* scene)
  { 
    BVH8* accel = new BVH8(Triangle8Type::type,scene);
    Accel::Intersectors intersectors= BVH8Triangle8Intersectors(accel);
    
    Builder* builder = NULL;
    if      (g_tri_builder == "default"     ) builder = BVH8Triangle8Builder(accel,scene,0);
    else if (g_tri_builder == "spatialsplit") builder = BVH8Triangle8Builder(accel,scene,MODE_HIGH_QUALITY);
    else if (g_tri_builder == "objectsplit" ) builder = BVH8Triangle8Builder(accel,scene,0);
    else THROW_RUNTIME_ERROR("unknown builder "+g_tri_builder+" for BVH8<Triangle8>");
    
    return new AccelInstance(accel,builder,intersectors);
  }

  Accel* BVH8::BVH8Triangle8ObjectSplit(Scene* scene)
  {
    BVH8* accel = new BVH8(Triangle8Type::type,scene);
    Accel::Intersectors intersectors= BVH8Triangle8Intersectors(accel);
    Builder* builder = BVH8Triangle8Builder(accel,scene,0);
    return new AccelInstance(accel,builder,intersectors);
  }

  Accel* BVH8::BVH8Triangle8SpatialSplit(Scene* scene)
  {
    BVH8* accel = new BVH8(Triangle8Type::type,scene);
    Accel::Intersectors intersectors= BVH8Triangle8Intersectors(accel);
    Builder* builder = BVH8Triangle8Builder(accel,scene,MODE_HIGH_QUALITY);
    return new AccelInstance(accel,builder,intersectors);
  }
}


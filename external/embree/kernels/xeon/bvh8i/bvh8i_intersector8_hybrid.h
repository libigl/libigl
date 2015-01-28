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

#ifndef __EMBREE_BVH8I_INTERSECTOR8_HYBRID_H__
#define __EMBREE_BVH8I_INTERSECTOR8_HYBRID_H__

#include "bvh8i.h"
#include "../common/stack_item.h"
#include "../common/ray8.h"

namespace embree
{
    
  namespace isa
  {
    /*! BVH8i Traverser. Packet traversal implementation for a Quad BVH. */
template<typename TriangleIntersector8>    
class BVH8iIntersector8Hybrid
    {

      /* shortcuts for frequently used types */
      typedef typename TriangleIntersector8::Primitive Triangle;
      typedef typename BVH4i::NodeRef NodeRef;
      typedef typename BVH8i::Node Node;
      typedef StackItemT<NodeRef> StackItem;
      static const size_t stackSizeSingle = 1+3*BVH4i::maxDepth;
      static const size_t stackSizeChunk = 4*BVH4i::maxDepth+1;

      static void intersect1(const BVH8i* bvh, NodeRef root, const size_t k, Ray8& ray, const avx3f &ray_org, const avx3f &ray_dir, const avx3f &ray_rdir, const avxf &ray_tnear, const avxf &ray_tfar, const avx3i& nearXYZ);
      static bool occluded1 (const BVH8i* bvh, NodeRef root, const size_t k, Ray8& ray, const avx3f &ray_org, const avx3f &ray_dir, const avx3f &ray_rdir, const avxf &ray_tnear, const avxf &ray_tfar, const avx3i& nearXYZ);

    public:
      static void intersect(avxb* valid, BVH8i* bvh, Ray8& ray);
      static void occluded (avxb* valid, BVH8i* bvh, Ray8& ray);
    };
  }
}

#endif
  

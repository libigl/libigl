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

#pragma once

#include "bvh8.h"
#include "common/ray4.h"
#include "common/stack_item.h"

namespace embree
{
  namespace isa 
  {
    /*! BVH8 Hybrid Packet traversal implementation. Switched between packet and single ray traversal. */
    template<typename PrimitiveIntersector>
      class BVH8Intersector4Hybrid 
    {
      /* shortcuts for frequently used types */
      typedef typename PrimitiveIntersector::Precalculations Precalculations;
      typedef typename PrimitiveIntersector::Primitive Primitive;
      typedef typename BVH8::NodeRef NodeRef;
      typedef typename BVH8::Node Node;
      static const size_t stackSizeSingle = 1+3*BVH8::maxDepth;
      static const size_t stackSizeChunk = 4*BVH8::maxDepth+1;

    public:
      static void intersect1(const BVH8* bvh, NodeRef root, size_t k, Precalculations& pre, Ray4& ray, 
			     const sse3f& ray_org, const sse3f& ray_dir, const sse3f& ray_rdir, const ssef& ray_tnear, const ssef& ray_tfar);
      static bool occluded1 (const BVH8* bvh, NodeRef root, size_t k, Precalculations& pre, Ray4& ray, 
			     const sse3f& ray_org, const sse3f& ray_dir, const sse3f& ray_rdir, const ssef& ray_tnear, const ssef& ray_tfar);

      static void intersect(sseb* valid, BVH8* bvh, Ray4& ray);
      static void occluded (sseb* valid, BVH8* bvh, Ray4& ray);
    };
  }
}

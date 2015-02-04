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

#include "bvh4.h"
#include "common/ray8.h"

namespace embree
{
  namespace isa
  {
    /*! BVH4 packet traversal implementation. */
    template<int types, bool robust, typename PrimitiveIntersector>
      class BVH4Intersector8Chunk
    {
      /* shortcuts for frequently used types */
      typedef typename PrimitiveIntersector::Precalculations Precalculations;
      typedef typename PrimitiveIntersector::Primitive Primitive;
      typedef typename BVH4::NodeRef NodeRef;
      typedef typename BVH4::Node Node;
      static const size_t stackSize = 4*BVH4::maxDepth+1;
      
    public:
      static void intersect(avxb* valid, BVH4* bvh, Ray8& ray);
      static void occluded (avxb* valid, BVH4* bvh, Ray8& ray);
    };
  }
}

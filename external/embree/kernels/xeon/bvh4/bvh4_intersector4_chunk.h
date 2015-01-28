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

#ifndef __EMBREE_BVH4_INTERSECTOR4_CHUNK_H__
#define __EMBREE_BVH4_INTERSECTOR4_CHUNK_H__

#include "bvh4.h"
#include "common/ray4.h"

namespace embree
{
  namespace isa 
  {
    /*! BVH4 packet traversal implementation. */
    template<typename PrimitiveIntersector>
      class BVH4Intersector4Chunk
    {
      /* shortcuts for frequently used types */
      typedef typename PrimitiveIntersector::Primitive Primitive;
      typedef typename BVH4::NodeRef NodeRef;
      typedef typename BVH4::Node Node;
      static const size_t stackSize = 4*BVH4::maxDepth+1;
      
    public:
      static void intersect(sseb* valid, BVH4* bvh, Ray4& ray);
      static void occluded (sseb* valid, BVH4* bvh, Ray4& ray);
    };
  }
}

#endif
  

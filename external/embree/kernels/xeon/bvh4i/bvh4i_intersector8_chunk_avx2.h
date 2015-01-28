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

#ifndef __EMBREE_BVH4I_INTERSECTOR8_CHUNK__AVX2_H__
#define __EMBREE_BVH4I_INTERSECTOR8_CHUNK__AVX2_H__

#include "bvh4i.h"
#include "../common/stack_item.h"
#include "../common/ray8.h"

namespace embree
{
  /*! BVH4i Traverser. Packet traversal implementation for a Quad BVH. */
  class BVH4iIntersector8ChunkAVX2
  {
  public:
    static void intersect(avxb* valid, BVH4i* bvh, Ray8& ray);
    static void occluded (avxb* valid, BVH4i* bvh, Ray8& ray);
  };
}

#endif

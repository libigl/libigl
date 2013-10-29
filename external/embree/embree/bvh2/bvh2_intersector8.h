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

#ifndef __EMBREE_BVH2_INTERSECTOR8_CHUNK_H__
#define __EMBREE_BVH2_INTERSECTOR8_CHUNK_H__

#include "bvh2.h"
#include "../include/intersector8.h"

namespace embree
{
  /*! BVH2 Traverser. Packet traversal implementation for a binary BVH. */
  template<typename TriangleIntersector>
    class BVH2Intersector8Chunk : public Intersector8
  {
    /* shortcuts for frequently used types */
    typedef typename TriangleIntersector::Triangle Triangle;
    typedef typename BVH2::NodeRef NodeRef;
    typedef typename BVH2::Node Node;

  public:

    BVH2Intersector8Chunk (const BVH2* bvh) 
      : Intersector8((intersectFunc)intersect,(occludedFunc)occluded), bvh(bvh) {}

    static Intersector8* create(const Accel* bvh) { 
      return new BVH2Intersector8Chunk((const BVH2*)bvh); 
    }

    static void   intersect(const BVH2Intersector8Chunk* This, Ray8& ray, const __m256 valid);
    static __m256 occluded (const BVH2Intersector8Chunk* This, Ray8& ray, const __m256 valid);

  private:
    const BVH2* bvh;
  };
}

#endif

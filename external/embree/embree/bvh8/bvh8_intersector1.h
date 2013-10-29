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

#ifndef __EMBREE_BVH8_INTERSECTOR1_H__
#define __EMBREE_BVH8_INTERSECTOR1_H__

#include "bvh8.h"
#include "../include/intersector1.h"
#include "../common/stack_item.h"

namespace embree
{
  /*! BVH8 Traverser. Single ray traversal implementation for a Quad BVH. */
  template<typename TriangleIntersector>
  class BVH8Intersector1 : public Intersector1
  {
    /* shortcuts for frequently used types */
    typedef typename TriangleIntersector::Triangle Triangle;
    typedef typename BVH8::NodeRef NodeRef;
    typedef typename BVH8::Node Node;
    typedef StackItemT<size_t> StackItem;

  public:
    BVH8Intersector1 (const BVH8* bvh) 
      : Intersector1((intersectFunc)intersect,(occludedFunc)occluded), bvh(bvh) {}

    static Intersector1* create(const Accel* bvh) { 
      return new BVH8Intersector1((const BVH8*)bvh); 
    }

    static void intersect(const BVH8Intersector1* This, Ray& ray);
    static bool occluded (const BVH8Intersector1* This, Ray& ray);

  private:
    const BVH8* bvh;
  };
}

#endif

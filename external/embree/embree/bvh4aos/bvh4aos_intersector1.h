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

#ifndef __EMBREE_BVH4AOS_INTERSECTOR1_H__
#define __EMBREE_BVH4AOS_INTERSECTOR1_H__

#include "bvh4aos.h"
#include "../include/intersector1.h"

namespace embree
{
  /*! BVH4AOS Traverser. Single ray traversal implementation for a Quad BVH. */
  template<typename TriangleIntersector>
  class BVH4AOSIntersector1 : public Intersector1
  {
    /* shortcuts for frequently used types */
    typedef typename TriangleIntersector::Triangle Triangle;
    typedef typename BVH4AOS::NodeRef NodeRef;
    typedef typename BVH4AOS::Node Node;
    
 public:
    BVH4AOSIntersector1 (const BVH4AOS* bvh) 
      : Intersector1((intersectFunc)intersect,(occludedFunc)occluded), bvh(bvh) {}

    static Intersector1* create(const Accel* bvh) { 
      return new BVH4AOSIntersector1((const BVH4AOS*)bvh); 
    }

    static void intersect(const BVH4AOSIntersector1* This, Ray& ray);
    static bool occluded (const BVH4AOSIntersector1* This, Ray& ray);

  private:
    const BVH4AOS* bvh;
  };
}

#endif

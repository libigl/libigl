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

#ifndef __EMBREE_BVH4AOS_TRIANGLE1_INTERSECTOR16_HYBRID_H__
#define __EMBREE_BVH4AOS_TRIANGLE1_INTERSECTOR16_HYBRID_H__

#include "bvh4aos.h"
#include "../common/ray16.h"
#include "../include/intersector16.h"
#include "../geometry/triangles.h"

namespace embree
{
  /*! BVH4AOS Traverser. Single ray traversal implementation for a Quad BVH. */
  class BVH4AOSTriangle1Intersector16Hybrid : public Intersector16
  {
    /* shortcuts for frequently used types */
    typedef Triangle1 Triangle;
    typedef typename BVH4AOS::NodeRef NodeRef;
    typedef typename BVH4AOS::Node Node;
    
  public:
    BVH4AOSTriangle1Intersector16Hybrid (const BVH4AOS* bvh) 
      : Intersector16((intersectFunc)intersect,(occludedFunc)occluded), bvh(bvh) {}

    static Intersector16* create(const Accel* bvh) { 
      return new BVH4AOSTriangle1Intersector16Hybrid((const BVH4AOS*)bvh); 
    }

    static void intersect (const BVH4AOSTriangle1Intersector16Hybrid* This, Ray16& ray, const __mmask valid);
    static __mmask occluded (const BVH4AOSTriangle1Intersector16Hybrid* This, Ray16& ray, const __mmask valid);

  private:
    static void intersect1(const BVH4AOS* bvh, 
                           const Node*     const nodes,
                           const Triangle1 *const triangles,
                           unsigned rootNode, long idx, 
                           const mic_f org_xyz, const mic_f dir_xyz, 
                           const mic_f rdir_xyz, const mic_f org_rdir_xyz, 
                           const mic_f min_dist_xyz, const mic_f max_dist_xyz, 
                           Ray16& ray);
    
    static bool occluded1(const BVH4AOS* bvh, 
                          const Node*     const nodes,
                          const Triangle1 *const triangles,
                          unsigned rootNode, long idx, 
                          const mic_f org_xyz, const mic_f dir_xyz, 
                          const mic_f rdir_xyz, const mic_f org_rdir_xyz, 
                          const mic_f min_dist_xyz, const mic_f max_dist_xyz);

  private:
    const BVH4AOS* bvh;
  };
}

#endif

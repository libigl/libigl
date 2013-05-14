// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#ifndef __EMBREE_BVH4_INTERSECTOR_H__
#define __EMBREE_BVH4_INTERSECTOR_H__

#include "bvh4.h"
#include "../common/intersector.h"

namespace embree
{
  /*! BVH4 Traverser. Single ray traversal implementation for a Quad BVH. */
  template<typename TriangleIntersector>
  class BVH4Intersector : public Intersector
  {
    /* shortcuts for frequently used types */
    typedef typename TriangleIntersector::Triangle Triangle;
    typedef typename BVH4::Base Base;
    typedef typename BVH4::Node Node;
    
  public:
    BVH4Intersector (const Ref<BVH4>& bvh) : bvh(bvh) {}
    void intersect(const Ray& ray, Hit& hit) const;
    bool occluded (const Ray& ray) const;

  private:
    Ref<BVH4> bvh;
  };
}

#endif

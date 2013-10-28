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

#ifndef __EMBREE_BVH_SORT_TEMPLATED_H__
#define __EMBREE_BVH_SORT_TEMPLATED_H__

#include "../common/default.h"
#include "../geometry/geometry.h"

namespace embree
{
  /* BVH child sorting by occlusion. */
  template<typename BVH>
    class BVHSortT 
  {
  public:
    typedef typename BVH::Node Node;
    typedef typename BVH::NodeRef NodeRef;

    static float sort(RTCGeometry* geom, BVH* bvh, NodeRef node)
    {
      /*! return if we found depth barrier */
      if (node.isBarrier())
        return 0;
      
      /*! sum up triangle areas */
      else if (node.isLeaf()) 
      {
        float A = 0.0f;
        size_t num; const char* tri = node.leaf(bvh->triPtr(),num);
        for (size_t i=0; i<num; i++)
          A += bvh->trity.area(tri+i*bvh->trity.bytes,geom);
        return A;
      }
      
      /*! sort node children based on approximated opacity */
      else 
      {
        /*! calculate half of node surface area */
        Node* n = node.node(bvh->nodePtr());
        
        /*! recurse into children first */
        float A = 0.0f;
        float opacity[BVH::N];
        for (size_t c=0; c<BVH::N; c++) {
          float dA = sort(geom,bvh,n->child(c));
          opacity[c] = dA*rcp(halfArea(n->bounds(c)));
          A += dA;
        }
      
        /*! sort children by opacity */
        for (size_t i=BVH::N-1; i>0; i--) {
          for (size_t j=0; j<i; j++) {
            if (opacity[j+0] > opacity[j+1]) {
              std::swap(opacity[j+0],opacity[j+1]);
              BVH::swap(n,j+0,n,j+1);
            }
          }
        }
        return A;
      }
    }
  };
}

#endif

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

#ifndef __EMBREE_BVH_ROTATE_TEMPLATED_H__
#define __EMBREE_BVH_ROTATE_TEMPLATED_H__

#include "../common/default.h"

namespace embree
{
  /* BVH Tree Rotations. */
  template<typename BVH>
    class BVHRotateT 
  {
  public:
    typedef typename BVH::Node Node;
    typedef typename BVH::NodeRef NodeRef;

    static size_t rotate(BVH* bvh, NodeRef parentRef, size_t depth, size_t maxDepth)
    {
      /*! nothing to rotate if we reached a leaf node. */
      if (depth >= maxDepth) return 0;
      if (parentRef.isBarrier()) return 0;
      if (parentRef.isLeaf()) return 0;
      Node* parent = parentRef.node(bvh->nodePtr());
      
      /*! rotate all children first */
      int max_cdepth = 0;
      int cdepth[BVH::N];
      for (size_t c=0; c<BVH::N; c++) {
        int d = (int)rotate(bvh,parent->child(c),depth+1,maxDepth);
        max_cdepth = max(max_cdepth,d);
        cdepth[c] = d;
      }
      
      /* compute bounds */
      float parentAreas[BVH::N]; BBox3f parentBounds[BVH::N];
      for (size_t i=0; i<BVH::N; i++) {
        const BBox3f bound = parent->bounds(i);
        parentBounds[i] = bound;
        parentAreas [i] = halfArea(bound);
      }
      
      /*! Find best rotation. We pick a first child a second child and a 
        subchild of it. We perform the best swap of the subchild and the 
        first child. */
      float bestArea = 0.0f;
      int bestChild1 = -1, bestChild2 = -1, bestChild2Child = -1;
      
      /* iterate over all first children */
      for (size_t c1=0; c1<BVH::N; c1++)
      {
        if ((depth+1)+cdepth[c1] > BVH::maxBuildDepth) continue;
        const BBox3f  child1Bound = parentBounds[c1];
        
        /* iterate over all second children */
        for (size_t c2=0; c2<BVH::N; c2++)
        {
          if (c1 == c2) continue;
          
          /*! ignore leaf nodes as we cannot descent into them */
          const NodeRef child2Ref = parent->child(c2);
          if (child2Ref.isBarrier()) continue;
          if (child2Ref.isLeaf()) continue;
          Node* child2 = child2Ref.node(bvh->nodePtr());
          
          /* compute child bounds */
          BBox3f child2Bounds[BVH::N];
          for (size_t i=0; i<BVH::N; i++)
            child2Bounds[i] = child2->bounds(i);
        
          /* iterate over all subchildren of the second child */
          for (size_t c2c=0; c2c<BVH::N; c2c++)
          {
            /* compute heuristic */
            BBox3f bounds = child1Bound;
            for (size_t i=0; i<BVH::N; i++) {
              if (i == c2c) continue;
              bounds = merge(bounds,child2Bounds[i]);
            }
            float area = halfArea(bounds)-parentAreas[c2];
            
            /* select best swap */
            if (area < bestArea) {
              bestArea = area;
              bestChild1 = (int)c1;
              bestChild2 = (int)c2;
              bestChild2Child = (int)c2c;
            }
          }
        }
      }
      
      /*! if we did not find a swap that improves the SAH then do nothing */
      if (bestChild1 == -1) return 1+max_cdepth;
      
      /*! perform the best found tree rotation */
      Node* child2 = parent->child(bestChild2).node(bvh->nodePtr());
      BVH::swap(parent,bestChild1,child2,bestChild2Child);
      parent->set(bestChild2,child2->bounds());
      
      /*! This returned depth is conservative as the child that was
       *  pulled up in the tree could have been on the critical path. */
      cdepth[bestChild1]++; // bestChild1 was pushed down one level
      return 1+max_cdepth; 
    }
  };
}

#endif

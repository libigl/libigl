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

#include "bvh4i.h"

namespace embree
{
  /* BVH4 Tree Rotations. */
  class BVH4iRotate
  {
  public:
    typedef BVH4i::Node Node;
    typedef BVH4i::NodeRef NodeRef;


  __forceinline static void mergedHalfArea( const BVH4i::Node *__restrict__ const parent,
					    const BVH4i::Node *__restrict__ const child2,
					    const size_t parentIndex,
					    float &cost,
					    int &pos) 
  { 
    const mic_f parent_min = parent->lowerXYZ(parentIndex);
    const mic_f parent_max = parent->upperXYZ(parentIndex);

    const mic_m lane0 = 0xf;
    const mic_m lane1 = 0xf0;
    const mic_m lane2 = 0xf00;
    const mic_m lane3 = 0xf000;

    const mic_f child2c0_min = select(lane0,parent_min,child2->lowerXYZ(0));
    const mic_f child2c0_max = select(lane0,parent_max,child2->upperXYZ(0));

    const mic_f child2c1_min = select(lane1,parent_min,child2->lowerXYZ(1));
    const mic_f child2c1_max = select(lane1,parent_max,child2->upperXYZ(1));

    const mic_f child2c2_min = select(lane2,parent_min,child2->lowerXYZ(2));
    const mic_f child2c2_max = select(lane2,parent_max,child2->upperXYZ(2));

    const mic_f child2c3_min = select(lane3,parent_min,child2->lowerXYZ(3));
    const mic_f child2c3_max = select(lane3,parent_max,child2->upperXYZ(3));

    const mic_f merged_min = min(min(child2c0_min,child2c1_min),min(child2c2_min,child2c3_min));
    const mic_f merged_max = max(max(child2c0_max,child2c1_max),max(child2c2_max,child2c3_max));

    const mic_f diag = merged_max - merged_min;
    const mic_f dx = swAAAA(diag);
    const mic_f dy = swBBBB(diag);
    const mic_f dz = swCCCC(diag);

    const mic_f half_area = dx*(dy+dz)+dy*dz; 

    const mic_f min_area = set_min_lanes(half_area);

    const mic_m m_min_area = eq(0x1111,min_area,half_area);
    const unsigned int i_min_area = bitscan(m_min_area);

    pos = i_min_area >> 2;
    cost = half_area[i_min_area];
  }


  public:
    static size_t rotate(BVH4i* bvh, NodeRef parentRef, size_t depth = 1,const bool onlyTopLevel = false);


    __forceinline static void rotateNode(BVH4i* bvh, NodeRef parentRef)
    {
      /*! nothing to rotate if we reached a leaf node. */
      if (parentRef.isLeaf()) return;
      assert(parentRef != BVH4i::invalidNode);

      Node* parent = parentRef.node(bvh->nodePtr());
      parent->prefetchNode<PFHINT_L1>();

      const size_t parentChildren = parent->numChildren();

      /* compute current area of all children */

      float childArea[4];

      while(1)
	{
	  const mic_f _childArea = parent->halfAreaBounds();
	  compactustore16f_low(0x1111,childArea,_childArea);

	  /*! Find best rotation. We pick a first child (child1) and a sub-child 
	    (child2child) of a different second child (child2), and swap child1 
	    and child2child. We perform the best such swap. */
	  float bestArea = 0.0f;
	  int bestChild1 = -1, bestChild2 = -1, bestChild2Child = -1;
	  for (size_t c2=0; c2<4; c2++)
	    {
	      /*! ignore leaf nodes as we cannot descent into them */
	      if (unlikely(parent->child(c2).isLeaf())) continue;
	      assert(parent->child(c2) != BVH4i::invalidNode);
	
	      Node* child2 = parent->child(c2).node(bvh->nodePtr());
	      child2->prefetchNode<PFHINT_L1>();

      
	      for (size_t i=0;i<parentChildren;i++)      
		{
		  if (unlikely(i == c2)) continue;
		  float area;
		  int pos;
		  mergedHalfArea(parent,child2,i,area,pos);

		  const float area_i = area - childArea[c2];

		  /*! accept a swap when it reduces cost and is not swapping a node with itself */
		  if (unlikely(area_i < bestArea)) {
		    bestArea = area_i;
		    bestChild1 = i;
		    bestChild2 = c2;
		    bestChild2Child = pos;
		  }
	  
		}
	    }    

	  /*! if we did not find a swap that improves the SAH then do nothing */
	  if (bestChild1 == -1) return;
    
	  /*! perform the best found tree rotation */
	  Node* child2 = parent->child(bestChild2).node(bvh->nodePtr());
	  assert( parent->child(bestChild1) != BVH4i::invalidNode);
	  assert( parent->child(bestChild2) != BVH4i::invalidNode);
	  assert( child2->child(bestChild2Child) != BVH4i::invalidNode);

	  BVH4i::swap(parent,bestChild1,child2,bestChild2Child);


	  parent->setBounds(bestChild2,child2);

	  BVH4i::compact(parent);
	  BVH4i::compact(child2);
	}
    }

  };
}

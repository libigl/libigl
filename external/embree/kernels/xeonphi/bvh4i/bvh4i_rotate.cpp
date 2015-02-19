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

#include "bvh4i_rotate.h"

#define DBG(x) 

#define FAST_ROTATE

namespace embree
{


  size_t BVH4iRotate::rotate(BVH4i* bvh, NodeRef parentRef, size_t depth, const bool onlyTopLevel)
  {
    DBG( DBG_PRINT(parentRef) );
    DBG( DBG_PRINT(depth) );

    assert(depth < BVH4i::maxBuildDepth);

    /*! nothing to rotate if we reached a leaf node. */
    if (parentRef.isLeaf()) return 0;
    assert(parentRef != BVH4i::invalidNode);

    Node* parent = parentRef.node(bvh->nodePtr());
    parent->prefetchNode<PFHINT_L1>();

    /*! rotate all children first */

    size_t cdepth[4];
    for (size_t c=0; c<4; c++)
      cdepth[c] = 0;

    const size_t parentChildren = parent->numChildren();

    for (size_t c=0; c<parentChildren; c++)
      {
	if (onlyTopLevel)
	  if (parent->upper[c].child == BVH4I_TOP_LEVEL_MARKER) continue; 

	cdepth[c] = rotate(bvh,parent->child(c),depth+1,onlyTopLevel);
      }

    /* compute current area of all children */
    float childArea[4];

#if defined(FAST_ROTATE)
    const mic_f _childArea = parent->halfAreaBounds();
    compactustore16f_low(0x1111,childArea,_childArea);
#else
    for (size_t i=0;i<4;i++)
      {
	childArea[i] = parent->halfAreaBounds( i );
	if (childArea[i] != _childArea[4*i])
	  FATAL("HERE");
	DBG( DBG_PRINT(i) );
	DBG( DBG_PRINT(childArea[i]) );

      }

    /*! get node bounds */
    BBox3fa child1_0,child1_1,child1_2,child1_3;
    child1_0 = parent->bounds( 0 );
    child1_1 = parent->bounds( 1 );
    child1_2 = parent->bounds( 2 );
    child1_3 = parent->bounds( 3 );
#endif

    /*! Find best rotation. We pick a first child (child1) and a sub-child 
      (child2child) of a different second child (child2), and swap child1 
      and child2child. We perform the best such swap. */
    float bestArea = 0;
    int bestChild1 = -1;
    int bestChild2 = -1;
    int bestChild2Child = -1;
    for (size_t c2=0; c2<4; c2++)
      {
	/*! ignore leaf nodes as we cannot descent into them */
	if (unlikely(parent->child(c2).isLeaf())) continue;
	assert(parent->child(c2) != BVH4i::invalidNode);

	Node* child2 = parent->child(c2).node(bvh->nodePtr());
	child2->prefetchNode<PFHINT_L1>();

#if !defined(FAST_ROTATE)
	/*! transpose child bounds */
	BBox3fa child2c0,child2c1,child2c2,child2c3;
	child2c0 = child2->bounds( 0 );
	child2c1 = child2->bounds( 1 );
	child2c2 = child2->bounds( 2 );
	child2c3 = child2->bounds( 3 );


	/*! put child1_0 at each child2 position */
	__aligned(16) float cost0[4];
	cost0[0] = halfArea(merge(child1_0,child2c1,child2c2,child2c3));
	cost0[1] = halfArea(merge(child2c0,child1_0,child2c2,child2c3));
	cost0[2] = halfArea(merge(child2c0,child2c1,child1_0,child2c3));
	cost0[3] = halfArea(merge(child2c0,child2c1,child2c2,child1_0));
	const float min0 = min(cost0[0],cost0[1],cost0[2],cost0[3]);
	int pos0 = (int)__bsf(mic_f(min0) == broadcast4to16f(cost0));
	assert(0 <= pos0 && pos0 < 4);

	float test_cost;
	int test_pos;

	mergedHalfArea(parent,child2,0,test_cost,test_pos);

	assert(min0 == test_cost);
	assert(pos0 == test_pos);

	/*! put child1_1 at each child2 position */
	__aligned(16) float cost1[4];
	cost1[0] = halfArea(merge(child1_1,child2c1,child2c2,child2c3));
	cost1[1] = halfArea(merge(child2c0,child1_1,child2c2,child2c3));
	cost1[2] = halfArea(merge(child2c0,child2c1,child1_1,child2c3));
	cost1[3] = halfArea(merge(child2c0,child2c1,child2c2,child1_1));

	const float min1 = min(cost1[0],cost1[1],cost1[2],cost1[3]);
	int pos1 = (int)__bsf(mic_f(min1) == broadcast4to16f(cost1));
	assert(0 <= pos1 && pos1 < 4);

	mergedHalfArea(parent,child2,1,test_cost,test_pos);
	assert(min1 == test_cost);
	assert(pos1 == test_pos);

	/*! put child1_2 at each child2 position */
	__aligned(16) float cost2[4];
	cost2[0] = halfArea(merge(child1_2,child2c1,child2c2,child2c3));
	cost2[1] = halfArea(merge(child2c0,child1_2,child2c2,child2c3));
	cost2[2] = halfArea(merge(child2c0,child2c1,child1_2,child2c3));
	cost2[3] = halfArea(merge(child2c0,child2c1,child2c2,child1_2));
	const float min2 = min(cost2[0],cost2[1],cost2[2],cost2[3]);
	int pos2 = (int)__bsf(mic_f(min2) == broadcast4to16f(cost2));
	assert(0 <= pos2 && pos2 < 4);

	mergedHalfArea(parent,child2,2,test_cost,test_pos);
	assert(min2 == test_cost);
	assert(pos2 == test_pos);

	/*! put child1_3 at each child2 position */
	__aligned(16) float cost3[4];
	cost3[0] = halfArea(merge(child1_3,child2c1,child2c2,child2c3));
	cost3[1] = halfArea(merge(child2c0,child1_3,child2c2,child2c3));
	cost3[2] = halfArea(merge(child2c0,child2c1,child1_3,child2c3));
	cost3[3] = halfArea(merge(child2c0,child2c1,child2c2,child1_3));
	const float min3 = min(cost3[0],cost3[1],cost3[2],cost3[3]);
	int pos3 = (int)__bsf(mic_f(min3) == broadcast4to16f(cost3));
	assert(0 <= pos3 && pos3 < 4);

	mergedHalfArea(parent,child2,3,test_cost,test_pos);
	assert(min3 == test_cost);
	assert(pos3 == test_pos);

	/*! find best other child */
	float area0123[4] = { min0,min1,min2,min3 };
	int pos[4] = { pos0,pos1,pos2,pos3 };

	for (size_t i=0;i<parentChildren;i++)
	  {	  
	    const float area_i = area0123[i] - childArea[c2];

	    DBG( DBG_PRINT(i) );
	    DBG( DBG_PRINT(area_i) );
	    DBG( DBG_PRINT(bestArea) );

	    if ((depth+1)+cdepth[i] > BVH4i::maxBuildDepth) continue;
	    if (i == c2) continue;

	    /*! accept a swap when it reduces cost and is not swapping a node with itself */
	    if (area_i < bestArea) {
	      bestArea = area_i;
	      bestChild1 = i;
	      bestChild2 = c2;
	      bestChild2Child = pos[i];

	      DBG( DBG_PRINT(bestChild1) );
	      DBG( DBG_PRINT(bestChild2) );
	    }
	  
	  }

#else
      
	for (size_t i=0;i<parentChildren;i++)
	  {	  
	    float area;
	    int pos;

	    mergedHalfArea(parent,child2,i,area,pos);

	    const float area_i = area - childArea[c2];

	    if ((depth+1)+cdepth[i] > BVH4i::maxBuildDepth) continue;
	    if (i == c2) continue;

	    /*! accept a swap when it reduces cost and is not swapping a node with itself */
	    if (area_i < bestArea) {
	      bestArea = area_i;
	      bestChild1 = i;
	      bestChild2 = c2;
	      bestChild2Child = pos;
	    }
	  
	  }
#endif

      }

    /*! if we did not find a swap that improves the SAH then do nothing */
    if (bestChild1 == -1) return 1+max(cdepth[0],cdepth[1],cdepth[2],cdepth[3]);
    
    /*! perform the best found tree rotation */
    Node* child2 = parent->child(bestChild2).node(bvh->nodePtr());
    assert( parent->child(bestChild1) != BVH4i::invalidNode);
    assert( parent->child(bestChild2) != BVH4i::invalidNode);
    assert( child2->child(bestChild2Child) != BVH4i::invalidNode);


#if DEBUG
    const float old_parent_sah = halfArea(parent->bounds());
    const float old_child2_sah = halfArea(child2->bounds());
#endif

    BVH4i::swap(parent,bestChild1,child2,bestChild2Child);


#if defined(FAST_ROTATE)
    parent->setBounds(bestChild2,child2);
#else
    parent->setBounds(bestChild2,child2->bounds());
#endif

    BVH4i::compact(parent);
    BVH4i::compact(child2);

#if DEBUG
    const float new_parent_sah = halfArea(parent->bounds());
    const float new_child2_sah = halfArea(child2->bounds());

    assert( new_parent_sah <= old_parent_sah);
    assert( new_child2_sah <= old_child2_sah);
    assert( parent->numChildren() >= 2 && parent->numChildren() <= 4);
    assert( child2->numChildren() >= 2 && child2->numChildren() <= 4);

#endif
    
    /*! This returned depth is conservative as the child that was
     *  pulled up in the tree could have been on the critical path. */
    cdepth[bestChild1]++; // bestChild1 was pushed down one level
    return 1+max(cdepth[0],cdepth[1],cdepth[2],cdepth[3]); 
  }
}

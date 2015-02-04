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

#include "bvh4_rotate.h"

namespace embree
{
  namespace isa 
  {
    /*! Computes half surface area of box. */
    __forceinline float halfArea3f(const BBox<ssef>& box) {
      const ssef d = box.size();
      const ssef a = d*shuffle<1,2,0,3>(d);
      return a[0]+a[1]+a[2];
    }
    
    size_t BVH4Rotate::rotate(BVH4* bvh, NodeRef parentRef, size_t depth)
    {
      /*! nothing to rotate if we reached a leaf node. */
      if (parentRef.isBarrier()) return 0;
      if (parentRef.isLeaf()) return 0;
      Node* parent = parentRef.node();
      
      /*! rotate all children first */
      ssei cdepth;
      for (size_t c=0; c<4; c++)
	cdepth[c] = (int)rotate(bvh,parent->child(c),depth+1);
      
      /* compute current area of all children */
      ssef sizeX = parent->upper_x-parent->lower_x;
      ssef sizeY = parent->upper_y-parent->lower_y;
      ssef sizeZ = parent->upper_z-parent->lower_z;
      ssef childArea = sizeX*(sizeY + sizeZ) + sizeY*sizeZ;
      
      /*! get node bounds */
      BBox<ssef> child1_0,child1_1,child1_2,child1_3;
      parent->bounds(child1_0,child1_1,child1_2,child1_3);
      
      /*! Find best rotation. We pick a first child (child1) and a sub-child 
	(child2child) of a different second child (child2), and swap child1 
	and child2child. We perform the best such swap. */
      float bestArea = 0;
      int bestChild1 = -1, bestChild2 = -1, bestChild2Child = -1;
      for (size_t c2=0; c2<4; c2++)
      {
	/*! ignore leaf nodes as we cannot descent into them */
	if (parent->child(c2).isBarrier()) continue;
	if (parent->child(c2).isLeaf()) continue;
	Node* child2 = parent->child(c2).node();
	
	/*! transpose child bounds */
	BBox<ssef> child2c0,child2c1,child2c2,child2c3;
	child2->bounds(child2c0,child2c1,child2c2,child2c3);
	
	/*! put child1_0 at each child2 position */
	float cost00 = halfArea3f(merge(child1_0,child2c1,child2c2,child2c3));
	float cost01 = halfArea3f(merge(child2c0,child1_0,child2c2,child2c3));
	float cost02 = halfArea3f(merge(child2c0,child2c1,child1_0,child2c3));
	float cost03 = halfArea3f(merge(child2c0,child2c1,child2c2,child1_0));
	ssef cost0 = ssef(cost00,cost01,cost02,cost03);
	ssef min0 = vreduce_min(cost0);
	int pos0 = (int)__bsf(movemask(min0 == cost0));
	
	/*! put child1_1 at each child2 position */
	float cost10 = halfArea3f(merge(child1_1,child2c1,child2c2,child2c3));
	float cost11 = halfArea3f(merge(child2c0,child1_1,child2c2,child2c3));
	float cost12 = halfArea3f(merge(child2c0,child2c1,child1_1,child2c3));
	float cost13 = halfArea3f(merge(child2c0,child2c1,child2c2,child1_1));
	ssef cost1 = ssef(cost10,cost11,cost12,cost13);
	ssef min1 = vreduce_min(cost1);
	int pos1 = (int)__bsf(movemask(min1 == cost1));
	
	/*! put child1_2 at each child2 position */
	float cost20 = halfArea3f(merge(child1_2,child2c1,child2c2,child2c3));
	float cost21 = halfArea3f(merge(child2c0,child1_2,child2c2,child2c3));
	float cost22 = halfArea3f(merge(child2c0,child2c1,child1_2,child2c3));
	float cost23 = halfArea3f(merge(child2c0,child2c1,child2c2,child1_2));
	ssef cost2 = ssef(cost20,cost21,cost22,cost23);
	ssef min2 = vreduce_min(cost2);
	int pos2 = (int)__bsf(movemask(min2 == cost2));
	
	/*! put child1_3 at each child2 position */
	float cost30 = halfArea3f(merge(child1_3,child2c1,child2c2,child2c3));
	float cost31 = halfArea3f(merge(child2c0,child1_3,child2c2,child2c3));
	float cost32 = halfArea3f(merge(child2c0,child2c1,child1_3,child2c3));
	float cost33 = halfArea3f(merge(child2c0,child2c1,child2c2,child1_3));
	ssef cost3 = ssef(cost30,cost31,cost32,cost33);
	ssef min3 = vreduce_min(cost3);
	int pos3 = (int)__bsf(movemask(min3 == cost3));
	
	/*! find best other child */
	ssef area0123 = ssef(extract<0>(min0),extract<0>(min1),extract<0>(min2),extract<0>(min3)) - ssef(childArea[c2]);
	int pos[4] = { pos0,pos1,pos2,pos3 };
	sseb valid = ssei(int(depth+1))+cdepth <= ssei(BVH4::maxBuildDepth); // only select swaps that fulfill depth constraints
	valid &= ssei(c2) != ssei(step);
	if (none(valid)) continue;
	size_t c1 = select_min(valid,area0123);
	float area = area0123[c1]; 
	assert(c1 != c2);
	
	/*! accept a swap when it reduces cost and is not swapping a node with itself */
	if (area < bestArea) {
	  bestArea = area;
	  bestChild1 = c1;
	  bestChild2 = c2;
	  bestChild2Child = pos[c1];
	}
      }
      
      /*! if we did not find a swap that improves the SAH then do nothing */
      if (bestChild1 == -1) return 1+reduce_max(cdepth);
      
      /*! perform the best found tree rotation */
      Node* child2 = parent->child(bestChild2).node();
      BVH4::swap(parent,bestChild1,child2,bestChild2Child);
      parent->set(bestChild2,child2->bounds());
      BVH4::compact(parent);
      BVH4::compact(child2);
      
      /*! This returned depth is conservative as the child that was
       *  pulled up in the tree could have been on the critical path. */
      cdepth[bestChild1]++; // bestChild1 was pushed down one level
      return 1+reduce_max(cdepth); 
    }
    
#if defined(__AVX__)
    
    struct RestructureNode
    {
      struct Partition3
      {
        __forceinline Partition3 () 
          : set0(0), set1(0), set2(0), sah(inf) {}
      public:
        int set0,set1,set2;
        float sah;
      };

      struct Item 
      {
        __forceinline Item () {}
        __forceinline Item (BVH4::NodeRef ref, const BBox3fa& b) 
          : ref(ref), bounds(b) { bounds.lower.a = ref.isLeaf() ? 0.0f : area(b); }

        __forceinline friend bool operator<(const Item& a, const Item& b) { 
          return a.bounds.lower.a < b.bounds.lower.a; 
        }
        
      public:
        BBox3fa bounds;
        BVH4::NodeRef ref;
      };

      RestructureNode (BVH4::Node* node2) 
      {
        for (size_t i=0; i< (1<<12); i++) 
          sa[i] = inf;

        numItems = 0;
        for (size_t c=0; c<4; c++) {
          if (node2->child(c) == BVH4::emptyNode) continue;
          items[numItems++] = Item(node2->child(c),node2->bounds(c));
        }
        if (numItems <= 1) return;
        std::sort(&items[0],&items[numItems]);

        BVH4::NodeRef ref0 = items[--numItems].ref;
        if (ref0.isLeaf()) return;
        BVH4::NodeRef ref1 = items[--numItems].ref;
        if (ref1.isLeaf()) return;
    
        BVH4::Node* node0 = ref0.node();
        for (size_t c=0; c<4; c++) {
          if (node0->child(c) == BVH4::emptyNode) continue;
          items[numItems++] = Item(node0->child(c),node0->bounds(c));
        }

        BVH4::Node* node1 = ref1.node();
        for (size_t c=0; c<4; c++) {
          if (node1->child(c) == BVH4::emptyNode) continue;
          items[numItems++] = Item(node1->child(c),node1->bounds(c));
        }

        Partition3 best;
        find012(0,0,0,0,best);

        int i0 = 0; node0->clear();
        int set0 = best.set0;
        int set1 = best.set1;
        int set2 = best.set2;
        while (set0) {
          const size_t j = __bscf(set0);
          node0->set(i0++,items[j].bounds,items[j].ref);
        }

        int i1 = 0; node1->clear();
        while (set1) {
          const size_t j = __bscf(set1);
          node1->set(i1++,items[j].bounds,items[j].ref);
        }

        int i2 = 0; node2->clear();
        while (set2) {
          const size_t j = __bscf(set2);
          node2->set(i2++,items[j].bounds,items[j].ref);
        }
        if (best.set0) {
          if (__popcnt(best.set0) == 1) {
            const size_t j = __bsf(best.set0);
            node2->set(i2++,items[j].bounds,items[j].ref);
          } else {
            node2->set(i2++,setBounds(best.set0),ref0);
          }
        }
        if (best.set1) {
          if (__popcnt(best.set1) == 1) {
            const size_t j = __bsf(best.set1);
            node2->set(i2++,items[j].bounds,items[j].ref);
          } else {
            node2->set(i2++,setBounds(best.set1),ref1);
          }
        }
      }

      void find012(int i, int set0, int set1, int set2, Partition3& best) 
      {
        if (i == numItems) return selectBest(set0,set1,set2,best);
        if (__popcnt(set0) == 4) return find12(i,set0,set1,set2,best);
        if (__popcnt(set1) == 4) return find02(i,set0,set1,set2,best);
        if (__popcnt(set2) == 2) return find01(i,set0,set1,set2,best);
        find012(i+1,set0|(1<<i),set1,set2,best);
        find012(i+1,set0,set1|(1<<i),set2,best);
        find012(i+1,set0,set1,set2|(1<<i),best);
      }

      void find01(int i, int set0, int set1, int set2, Partition3& best) 
      {
        if (i == numItems) return selectBest(set0,set1,set2,best);
        if (__popcnt(set0) == 4) return find1(i,set0,set1,set2,best);
        if (__popcnt(set1) == 4) return find0(i,set0,set1,set2,best);
        find01(i+1,set0|(1<<i),set1,set2,best);
        find01(i+1,set0,set1|(1<<i),set2,best);
      }

      void find02(int i, int set0, int set1, int set2, Partition3& best) 
      {
        if (i == numItems) return selectBest(set0,set1,set2,best);
        if (__popcnt(set0) == 4) return find2(i,set0,set1,set2,best);
        if (__popcnt(set2) == 2) return find0(i,set0,set1,set2,best);
        find02(i+1,set0|(1<<i),set1,set2,best);
        find02(i+1,set0,set1,set2|(1<<i),best);
      }

      void find12(int i, int set0, int set1, int set2, Partition3& best) 
      {
        if (i == numItems) return selectBest(set0,set1,set2,best);
        if (__popcnt(set1) == 4) return find2(i,set0,set1,set2,best);
        if (__popcnt(set2) == 2) return find1(i,set0,set1,set2,best);
        find12(i+1,set0,set1|(1<<i),set2,best);
        find12(i+1,set0,set1,set2|(1<<i),best);
      }
      
      void find0(int i, int set0, int set1, int set2, Partition3& best) 
      {
        if (i == numItems) return selectBest(set0,set1,set2,best);
        find0(i+1,set0|(1<<i),set1,set2,best);
      }

      void find1(int i, int set0, int set1, int set2, Partition3& best) 
      {
        if (i == numItems) return selectBest(set0,set1,set2,best);
        find1(i+1,set0,set1|(1<<i),set2,best);
      }
      
      void find2(int i, int set0, int set1, int set2, Partition3& best) 
      {
        if (i == numItems) return selectBest(set0,set1,set2,best);
        find2(i+1,set0,set1,set2|(1<<i),best);
      }

      void selectBest(int set0, int set1, int set2, Partition3& best)
      {
        if (sa[set0] == float(inf)) sa[set0] = setArea(set0);
        if (sa[set1] == float(inf)) sa[set1] = setArea(set1);
        const float sah = sa[set0] + sa[set1];
        if (best.sah < sah) return;
        best.sah  = sah;
        best.set0 = set0;
        best.set1 = set1;
        best.set2 = set2;
      }
      
      float setArea(size_t set) 
      {
        if (set == 0) 
          return 0.0f;
        
        BBox3fa bounds = empty;
        while (set) {
          bounds.extend(items[__bscf(set)].bounds);
        }
        return area(bounds);
      }

      BBox3fa setBounds(size_t set) 
      {
        if (set == 0) 
          return empty;
        
        BBox3fa bounds = empty;
        while (set) {
          bounds.extend(items[__bscf(set)].bounds);
        }
        return bounds;
      }
      
    private:
      float sa[1<<12];
      Item items[12];
      size_t numItems;
    };

    void BVH4Rotate::restructure(NodeRef ref, size_t depth)
    {
      if (ref.isLeaf()) return;
      Node* node = ref.node();
      RestructureNode temp(node);
      for (size_t c=0; c<4; c++) {
        if (node->child(c) == BVH4::emptyNode) continue;
        restructure(node->child(c),depth+1);
      }
    }
#endif
  }
}

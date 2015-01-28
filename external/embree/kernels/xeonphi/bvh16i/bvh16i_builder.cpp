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

#include "bvh16i/bvh16i_builder.h"

#define ENABLE_EXTENDED_LEAVES

namespace embree
{
#define DBG(x) 

  static mic_i bvh16i_node_dist = 0;
  static mic_i bvh4i_node_dist = 0;

  Builder* BVH16iBuilder::create (void* accel, BuildSource* source, void* geometry, size_t mode ) 
  { 
    Builder* builder = new BVH16iBuilder((BVH16i*)accel,source,geometry);
    return builder;
  }

  void BVH16iBuilder::printBuilderName()
  {
    std::cout << "building BVH16i with binned SAH builder (MIC) ... " << std::endl;    
  }

  void BVH16iBuilder::convertQBVHLayout(const size_t threadIndex, const size_t threadCount)
  {
    DBG(PING);
    DBG_PRINT(numPrimitives);

    countLeaves(0);

    //BVH16i::Node *bvh16 = (BVH16i::Node*)this->prims;
    const size_t sizeBVH4i = atomicID * sizeof(BVHNode);
    BVH16i::Node *bvh16 = (BVH16i::Node*)(os_malloc(sizeof(BVH16i::Node)*numPrimitives));

    
    bvh16[0].reset();
    size_t index16 = 1;
    bvh16i_node_dist = 0;
    bvh4i_node_dist = 0;
    convertBVH4iToBVH16i(node,
			 node[0].lower.a,
			 node[0].upper.a,
			 bvh16,
			 index16,
			 (unsigned int&)bvh16[0].child[0]);
        
    DBG_PRINT(numPrimitives * sizeof(BVHNode) / sizeof(BVH16i::Node));

    const size_t sizeBVH16i = index16 * sizeof(BVH16i::Node);

    DBG_PRINT(index16);
    DBG_PRINT(sizeBVH16i);
    DBG_PRINT((float)sizeBVH4i / sizeBVH16i);
 
    /* bvh16i node util */
    {
      unsigned int total = 0;
      float util = 0.0f;
      for (size_t i=0;i<16;i++) {
	util += (float)(i+1) * bvh16i_node_dist[i];
	total += bvh16i_node_dist[i];
      }
      DBG_PRINT(total);
      std::cout << "bvh16i node util dist: ";
      DBG_PRINT(bvh16i_node_dist);
      float sum = 0;
      for (size_t i=0;i<16;i++) 
	{
	  sum += (float)bvh16i_node_dist[i] * 100.0f / total;
	  std::cout << i+1 << "[" << (float)bvh16i_node_dist[i] * 100.0f / total << "%, sum " << sum << "%] ";
	}
      std::cout << std::endl;
      DBG_PRINT(100.0f * util / (16.0f * total));
      std::cout << std::endl;
    }

    /* bvh4i node util */
    {
      unsigned int total = 0;
      float util = 0.0f;
      for (size_t i=0;i<16;i++) {
	util += (float)(i+1) * bvh4i_node_dist[i];
	total += bvh4i_node_dist[i];
      }
      DBG_PRINT(total);
      std::cout << "bvh4i node util dist: ";
      DBG_PRINT(bvh4i_node_dist);
      float sum = 0;
      for (size_t i=0;i<16;i++) 
	{
	  sum += (float)bvh4i_node_dist[i] * 100.0f / total;
	  std::cout << i+1 << "[" << (float)bvh4i_node_dist[i] * 100.0f / total << "%, sum " << sum << "%] ";
	}
      std::cout << std::endl;
      DBG_PRINT(100.0f * util / (16.0f * total));
      std::cout << std::endl;
    }


    ((BVH16i*)bvh)->node16 = bvh16;

    // for bvh4i
    LockStepTaskScheduler::dispatchTask( task_convertToSOALayout, this, threadIndex, threadCount );    
  }

  void BVH16iBuilder::countLeaves(const size_t index)
  {
    BVHNode &entry = node[index];

    if (entry.isLeaf())
      {
	entry.upper.a = 1;
      }
    else
      {
	const unsigned int childID = entry.firstChildID();
	const unsigned int children = entry.items();

	unsigned int leaves = 0;
	for (unsigned int i=0;i<children;i++) 
	  {
	    countLeaves(childID+i);
	    leaves += node[childID+i].upper.a;
	  }      
	entry.upper.a = leaves;
      }
  }

  void BVH16iBuilder::getLeaves(unsigned int bvh4_ext_min, 
				unsigned int node_index, 
				BVHNode *leaves, 
				unsigned int &numLeaves)
  {

    if (bvhLeaf(bvh4_ext_min))
      {
	leaves[numLeaves++] = node[node_index];
      }
    else
      {
        const size_t childID = bvhChildID(bvh4_ext_min);
        const size_t children = bvhItems(bvh4_ext_min);

	for (unsigned int i=0;i<children;i++) 
	  {
	    BVHNode &entry = node[childID+i];
	    getLeaves(entry.lower.a,childID+i,leaves,numLeaves);
	  }
      }
  }



  void BVH16iBuilder::convertBVH4iToBVH16i(const BVHNode *const bvh4,
					   const unsigned int bvh4_ext_min, 
					   const unsigned int numLeavesInSubTree, 
					   BVH16i::Node *const bvh16,
					   size_t &index16,
					   unsigned int &parent_offset)
    {
      size_t bvh16_used_slots = 0;
      const size_t bvh16_node_index = index16++;
      
      DBG(DBG_PRINT(bvh16_node_index));
      
#if defined(ENABLE_EXTENDED_LEAVES)
      if (numLeavesInSubTree <= 4)
       	{
       	  BVHNode *leaves = (BVHNode*)&bvh16[bvh16_node_index];
	  const mic_f init_node = load16f((float*)BVH4i::initQBVHNode);
	  store16f_ngo((float*)&leaves[0],init_node);
	  store16f_ngo((float*)&leaves[2],init_node);

	  unsigned int numLeaves = 0;
       	  getLeaves(bvh4_ext_min,(unsigned int)-1,leaves,numLeaves);

	  bvh4i_node_dist[numLeaves-1]++;

	  if (numLeavesInSubTree != numLeaves) FATAL("HERE");

	  convertToBVH4Layout(leaves);

	  parent_offset  = (unsigned int)(sizeof(BVH16i::Node) * bvh16_node_index);
	  parent_offset |= BVH16I_EXTENDED_LEAF_MASK;
	  
	  return;
	}
#endif		    


      bvh16[bvh16_node_index].reset();
      
      {

        const size_t childID = bvhChildID(bvh4_ext_min);
        const size_t children = bvhItems(bvh4_ext_min);

        DBG(
          DBG_PRINT(childID);
          DBG_PRINT(children);
          );
        
        for (size_t i=0;i<children;i++) 
	{
	  DBG(std::cout << "Putting " << bvh4[childID+i] << " in slot " << bvh16_used_slots << std::endl);
	  bvh16[bvh16_node_index].set(bvh16_used_slots++,bvh4[childID+i]);
	}
      }
      
      DBG(DBG_PRINT(bvh16[bvh16_node_index]));
      
      while(bvh16_used_slots < 16)
      {
	DBG(
          std::cout << std::endl << std::flush;
          DBG_PRINT(bvh16_used_slots);
          DBG_PRINT(bvh16[bvh16_node_index]);
          );
	mic_f node_area = bvh16[bvh16_node_index].area();
	DBG(DBG_PRINT(node_area));
        
        
	ssize_t max_index = -1;
	float max_area = 0.0f;
	const unsigned int free_slots = 16 - bvh16_used_slots;
        
	ssize_t max_index_small = -1;
	ssize_t min_children_small = 16;
	float max_area_small = 0.0f;
        
	for (size_t i=0;i<bvh16_used_slots;i++)
        {
          if (bvhLeaf(bvh16[bvh16_node_index].child[i])) continue;
          if ((bvhItems(bvh16[bvh16_node_index].child[i]) + bvh16_used_slots - 1) <= 16 && 
              node_area[i] > max_area)
          {      
            //if (bvh16_used_slots >=8)
            if (bvh16[bvh16_node_index].data[i] >= 8 && bvh16[bvh16_node_index].data[i] <= 16) continue;	      
            
            max_index = i;
            max_area = node_area[i];
            
            DBG(
              DBG_PRINT(i);
              DBG_PRINT(max_index);
              DBG_PRINT(max_area);
              );
            
          }
	  
          if (bvh16[bvh16_node_index].data[i] <= free_slots && 
              bvh16[bvh16_node_index].data[i] < min_children_small)// &&
            //node_area[i] > max_area_small)
          {
            min_children_small = bvh16[bvh16_node_index].data[i];
            max_index_small = i;
            max_area_small = node_area[i];
            
            DBG(
              DBG_PRINT(i);
              DBG_PRINT(max_index_small);
              DBG_PRINT(max_area_small);
              );
            
          }	  
          
        }
        
	if (max_index == -1) 
        {
          break;
        }
        
	if (max_index_small != -1)
	  max_index = max_index_small;
        
        
	DBG(std::cout << "1" << std::endl << std::flush);
	const unsigned int parent_index = bvh16[bvh16_node_index].child[max_index];
	DBG(DBG_PRINT(parent_index));
        
	
	bvh16[bvh16_node_index].shift(max_index);
	bvh16_used_slots--;
        
	assert(!bvhLeaf(parent_index));
	assert( bvhItems(parent_index) + bvh16_used_slots <= 16);
        
	DBG(std::cout << "2" << std::endl << std::flush);
        
	const size_t childID = bvhChildID(parent_index);
	const size_t children = bvhItems(parent_index);
	for (size_t i=0;i<children;i++) 
        {
          DBG(std::cout << "Putting node " << childID+i << " -> " << bvh4[childID+i] << " in slot " << bvh16_used_slots << std::endl);
          bvh16[bvh16_node_index].set(bvh16_used_slots++,bvh4[childID+i]);
          
        }
	if (bvh16_used_slots > 16) FATAL("HERE");
        
	DBG(std::cout << "3" << std::endl << std::flush);
	assert(bvh16_used_slots <= 16);
      }
      
      
      DBG(
        std::cout << "FINAL: " << std::endl << std::flush;
        DBG_PRINT(bvh16_used_slots);
        DBG_PRINT(bvh16[bvh16_node_index]);
        );
      
      DBG(DBG_PRINT(bvh16[bvh16_node_index].child));
      
      
      DBG(DBG_PRINT(bvh16_used_slots));
      
      parent_offset = (unsigned int)(sizeof(BVH16i::Node) * bvh16_node_index);
      
      bvh16i_node_dist[bvh16_used_slots-1]++;
      
      BVH16i::Node &b16 = bvh16[bvh16_node_index];


      DBG(DBG_PRINT(b16));
      for (size_t i=0;i<bvh16_used_slots;i++)
        if (!bvhLeaf(b16.child[i]))
	{
	  DBG(std::cout << "RECURSE FOR " << b16.child[i] << " " << b16.data[i] << std::endl << std::flush);

	  //DBG_PRINT(b16.data[i]); // TRY: EITHER BVH4i leaf/leaves or triangle16 struct 

	  convertBVH4iToBVH16i(bvh4,
			       b16.child[i],
			       b16.data[i],
			       bvh16,			      
			       index16,
			       (unsigned int&)b16.child[i]);
	}
        else
	{
	  b16.child[i] = (b16.child[i] ^ BVH_LEAF_MASK) | QBVH_LEAF_MASK;
	}
      
      DBG(std::cout << "DONE" << std::endl << std::flush);
      
    }


}

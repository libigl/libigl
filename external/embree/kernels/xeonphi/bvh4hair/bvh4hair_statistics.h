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

#include "bvh4hair.h"

namespace embree
{
  template<typename NodeType>
  class BVH4HairStatistics 
  {
    typedef NodeType Node;
    typedef BVH4Hair::NodeRef NodeRef;

  public:

    /* Constructor gathers statistics. */
    BVH4HairStatistics (BVH4Hair* bvh);

    /*! Convert statistics into a string */
    std::string str();

  private:
    void statistics(NodeRef node, BBox3fa bounds, size_t& depth);

  private:
    BVH4Hair* bvh;
    float bvhSAH;                      //!< SAH cost of the BVH4.
    float leafSAH;                     //!< SAH cost of the BVH4.
    size_t numNodes;                   //!< Number of internal nodes.
    size_t numValidBoxes;              //!< Number of valid boxes per node.
    size_t numAlignedNodes;            //!< Number of internal aligned nodes.
    size_t numUnalignedNodes;          //!< Number of internal aligned nodes.
    size_t numLeaves;                  //!< Number of leaf nodes.
    size_t numPrimBlocks;              //!< Number of primitive blocks.
    size_t numPrims;                   //!< Number of primitives.
    size_t depth;                      //!< Depth of the tree.
  };


  template<typename NodeType>
  BVH4HairStatistics<NodeType>::BVH4HairStatistics (BVH4Hair* bvh) : bvh(bvh)
  {
    numNodes = numLeaves = numPrimBlocks = numAlignedNodes = numUnalignedNodes = numPrims = depth = 0;
    numValidBoxes = 0;
    bvhSAH = leafSAH = 0.0f;
    if (bvh->root != BVH4Hair::invalidNode)
      statistics(bvh->root,bvh->bounds,depth);
    bvhSAH /= area(bvh->bounds);
    leafSAH /= area(bvh->bounds);
    assert(depth <= BVH4Hair::maxDepth);
  }

  template<typename NodeType>
  std::string BVH4HairStatistics<NodeType>::str()  
  {
    std::ostringstream stream;
    size_t bytesNodes = bvh->size_node; 
    size_t bytesPrims  = bvh->size_accel;
    size_t bytesTotal = bytesNodes+bytesPrims;
    size_t bytesTotalAllocated = bvh->bytes();
    stream.setf(std::ios::fixed, std::ios::floatfield);
    stream.precision(4);
    stream << "  sah = " << bvhSAH << ", leafSAH = " << leafSAH;
    stream.setf(std::ios::fixed, std::ios::floatfield);
    stream.precision(1);
    stream << ", depth = " << depth;
    stream << ", size = " << bytesTotalAllocated/1E6 << " MB " << std::endl;
    stream.precision(1);
    stream << "  nodes = "  << numNodes << " "
           << "(" << double(sizeof(NodeType)*numNodes)/1E6 << " MB used of " << bytesNodes/1E6  << " MB pre-allocated) "
           << "(" << 100.0*double(bytesNodes)/double(bytesTotal) << "% of total " << bytesTotalAllocated/1E6 << " MB) "
           << "(" << 100.0*double(sizeof(NodeType)*numNodes)/double(bytesNodes) << "% of allocated nodes used)" 
           << std::endl;
    stream << "  AABB nodes " << numAlignedNodes 
	   << " (" << 100.0*double(numAlignedNodes)/double(numNodes) << "% of total) "
	   << " / OBB nodes " << numUnalignedNodes
	   << " (" << 100.0*double(numUnalignedNodes)/double(numNodes) << "% of total) " << std::endl;
    stream << "  leaves = " << numLeaves << " "
           << "(" << bytesPrims/1E6  << " MB) "
           << "(" << 100.0*double(bytesPrims)/double(bytesTotal) << "% of total " << bytesTotalAllocated/1E6 << " MB) "
           << std::endl;
    stream << "  node utilization " << 100.0*double(numValidBoxes)/double(numNodes*4.0) << "%" << std::endl;
    stream << "  leaf utilization " << 100.0*double(numPrims)/double(numPrimBlocks*2) << "%" << std::endl;

    return stream.str();
  }

  template<typename NodeType>
  void BVH4HairStatistics<NodeType>::statistics(NodeRef node, BBox3fa bounds, size_t& depth)
  {
    float A = bounds.empty() ? 0.0f : area(bounds);
    
    if (node.isNode())
    {
      if (node & BVH4Hair::alignednode_mask) 
	numAlignedNodes++;
      else
	numUnalignedNodes++;
      node = (unsigned int)node & (~BVH4Hair::alignednode_mask);
      numNodes++;
      depth = 0;
      size_t cdepth = 0;
      Node* n = (Node*)node.node((Node*)bvh->nodePtr());

      bvhSAH += A*BVH4Hair::travCost;

      for (size_t i=0; i<BVH4Hair::N; i++) {
	if (n->child(i) == BVH4Hair::invalidNode) {continue; }
	numValidBoxes++;
        statistics(n->child(i),n->bounds(i),cdepth); 
        depth=max(depth,cdepth);
      }
      depth++;
      return;
    }
    else
    {
      depth = 0;
      unsigned int prims; const char* tri = node.leaf(bvh->primitivesPtr(),prims);
      if (!prims) return;
      
      numLeaves++;
      numPrimBlocks  += 1;
      numPrims       += prims;
      float sah = A * bvh->primTy.intCost * 1;
      bvhSAH += sah;
      leafSAH += sah;
    }
  }

}

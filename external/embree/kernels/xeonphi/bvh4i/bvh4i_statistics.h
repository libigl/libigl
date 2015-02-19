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
  template<typename NodeType>
  class BVH4iStatistics 
  {
    typedef NodeType Node;
    typedef BVH4i::NodeRef NodeRef;

  public:

    /* Constructor gathers statistics. */
    BVH4iStatistics (BVH4i* bvh);

    /*! Convert statistics into a string */
    std::string str();

  private:
    void statistics(NodeRef node, BBox3fa bounds, size_t& depth);

  private:
    BVH4i* bvh;
    float bvhSAH;                      //!< SAH cost of the BVH4.
    float leafSAH;                      //!< SAH cost of the BVH4.
    size_t numNodes;                   //!< Number of internal nodes.
    size_t numValidBoxes;              //!< Number of valid boxes per node.
    size_t numLeaves;                  //!< Number of leaf nodes.
    size_t numPrimBlocks;              //!< Number of primitive blocks.
    size_t numPrimBlocks4;             //!< Number of primitive blocks, assuming block size of 4
    size_t numPrims;                   //!< Number of primitives.
    size_t depth;                      //!< Depth of the tree.
  };


  template<typename NodeType>
  BVH4iStatistics<NodeType>::BVH4iStatistics (BVH4i* bvh) : bvh(bvh)
  {
    numNodes = numLeaves = numPrimBlocks = numPrimBlocks4 = numPrims = depth = 0;
    numValidBoxes = 0;
    bvhSAH = leafSAH = 0.0f;
    if (bvh->root != BVH4i::invalidNode)
      statistics(bvh->root,bvh->bounds,depth);
    bvhSAH /= area(bvh->bounds);
    leafSAH /= area(bvh->bounds);
    assert(depth <= BVH4i::maxDepth);
  }

  template<typename NodeType>
  std::string BVH4iStatistics<NodeType>::str()  
  {
    std::ostringstream stream;
    size_t bytesNodes = bvh->size_node; // numNodes*sizeof(Node);
    size_t bytesTris  = bvh->size_accel;
    size_t bytesTotal = bytesNodes+bytesTris;//+bytesVertices;
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
    stream << "  leaves = " << numLeaves << " "
           << "(" << bytesTris/1E6  << " MB) "
           << "(" << 100.0*double(bytesTris)/double(bytesTotal) << "% of total " << bytesTotalAllocated/1E6 << " MB) "
           << std::endl;
    stream << "  node utilization " << 100.0*double(numValidBoxes)/double(numNodes*4) << "%" << std::endl;
    stream << "  leaf utilization " << 100.0*double(numPrims)/double(numPrimBlocks*4) << "%" << std::endl;
    return stream.str();
  }

  template<typename NodeType>
  void BVH4iStatistics<NodeType>::statistics(NodeRef node, BBox3fa bounds, size_t& depth)
  {
    float A = bounds.empty() ? 0.0f : area(bounds);
    
    if (!isfinite(A))
      {
	DBG_PRINT(node);
	DBG_PRINT(bounds);
	DBG_PRINT(depth);
	FATAL("error in sah");
      }

    if (node.isNode())
    {
      numNodes++;
      depth = 0;
      size_t cdepth = 0;
      Node* n = (Node*)node.node((Node*)bvh->nodePtr());

      bvhSAH += A*BVH4i::travCost;

      for (size_t i=0; i<BVH4i::N; i++) {
	if (n->child(i) == BVH4i::invalidNode) { break; }
	numValidBoxes++;

	BBox3fa b = n->bounds(i);
	if (!(isfinite(b.lower.x) && isfinite(b.lower.y) && isfinite(b.lower.z)))
	  FATAL("lower");

	if (!(isfinite(b.upper.x) && isfinite(b.upper.y) && isfinite(b.upper.z)))
	  FATAL("upper");

        statistics(n->child(i),n->bounds(i),cdepth); 
        depth=max(depth,cdepth);
      }
      depth++;
      return;
    }
    else
    {
      depth = 0;
      unsigned int prims; const char* tri = node.leaf(bvh->triPtr(),prims);
      assert(prims > 0);
      if (!prims) return;
      
      numLeaves++;
      numPrimBlocks  += 1;
      numPrims       += prims;
      numPrimBlocks4 += (prims+3)/4;
      float sah = A * bvh->primTy.intCost * 1;
      bvhSAH += sah;
      leafSAH += sah;
    }
  }

}

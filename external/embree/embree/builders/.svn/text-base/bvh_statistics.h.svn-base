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

#ifndef __EMBREE_BVH_STATISTICS_TEMPLATED_H__
#define __EMBREE_BVH_STATISTICS_TEMPLATED_H__

#include "../common/default.h"

namespace embree
{
  /* BVH Tree Statistics. */
  template<typename BVH>
    class BVHStatisticsT 
  {
  public:
    typedef typename BVH::Node Node;
    typedef typename BVH::NodeRef NodeRef;

    BVHStatisticsT (BVH* bvh) : bvh(bvh)
    {
      /* calculate statistics */
      numNodes = numLeaves = numPrimBlocks = numPrims = depth = 0;
      bvhSAH = statistics(bvh->root,0.0f,depth);
      assert(depth <= BVH::maxDepth);
    }

    /*! Print statistics of the BVH. */
    void print()  
    {
      /* output statistics */
      std::ostringstream stream;
      size_t bytesNodes = numNodes     *sizeof(Node);
      size_t bytesTris  = numPrimBlocks*bvh->trity.bytes;
      size_t numVertices = bvh->numVertices;
      size_t bytesVertices = numVertices*sizeof(Vec3fa); 
      size_t bytesTotal = bytesNodes+bytesTris+bytesVertices;
      stream.setf(std::ios::scientific, std::ios::floatfield);
      stream.precision(2);
      stream << "sah = " << bvhSAH << std::endl;
      stream.setf(std::ios::fixed, std::ios::floatfield);
      stream.precision(1);
      stream << "depth = " << depth << std::endl;
      stream << "size = " << bytesTotal/1E6 << " MB" << std::endl;
      stream.precision(1);
      stream << "nodes = "  << numNodes << " "
             << "(" << bytesNodes/1E6  << " MB) "
             << "(" << 100.0*double(bytesNodes)/double(bytesTotal) << "% of total) "
             << "(" << 100.0*(numNodes-1+numLeaves)/(BVH::N*numNodes) << "% used)" 
             << std::endl;
      stream << "leaves = " << numLeaves << " "
             << "(" << bytesTris/1E6  << " MB) "
             << "(" << 100.0*double(bytesTris)/double(bytesTotal) << "% of total) "
             << "(" << 100.0*numPrims/(bvh->trity.blockSize*numPrimBlocks) << "% used)" 
             << std::endl;
      stream << "vertices = " << numVertices << " "
             << "(" << bytesVertices/1E6 << " MB) " 
             << "(" << 100.0*double(bytesVertices)/double(bytesTotal) << "% of total) "
             << "(" << 100.0*12.0f/float(sizeof(Vec3fa)) << "% used)" 
             << std::endl;
      std::cout << stream.str();
    }

  private:

    float statistics(NodeRef node, float ap, size_t& depth)
    {
      if (node.isNode())
      {
        numNodes++;
        depth = 0;
        size_t cdepth = 0;
        const Node* n = node.node(bvh->nodePtr());
        float sah = 0.0f;
        for (size_t i=0; i<BVH::N; i++) {
          sah += statistics(n->child(i),area(n->bounds(i)),cdepth); 
          depth=max(depth,cdepth);
        }
        depth++;
        return ap*BVH::travCost + sah;
      }
      else
      {
        depth = 0;
        size_t num; const char* tri = node.leaf(bvh->triPtr(),num);
        if (!num) return 0.0f;
        
        numLeaves++;
        numPrimBlocks += num;
        for (size_t i=0; i<num; i++)
          numPrims += bvh->trity.size(tri+i*bvh->trity.bytes);
        return bvh->trity.intCost * ap * num;
      }
    }

    BVH* bvh;
    float bvhSAH;                      //!< SAH cost of the BVH.
    size_t numNodes;                   //!< Number of internal nodes.
    size_t numLeaves;                  //!< Number of leaf nodes.
    size_t numPrimBlocks;              //!< Number of primitive blocks.
    size_t numPrims;                   //!< Number of primitives.
    size_t depth;                      //!< Depth of the tree.
  };
}

#endif

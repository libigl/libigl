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

#include "bvh4i_statistics.h"
#include "geometry/triangle1.h"

namespace embree
{
  BVH4iStatistics::BVH4iStatistics (BVH4i* bvh) : bvh(bvh)
  {
    numNodes = numLeaves = numPrimBlocks = numPrimBlocks4 = numPrims = depth = 0;
    bvhSAH = leafSAH = 0.0f;
    statistics(bvh->root,bvh->bounds,depth);
    bvhSAH /= area(bvh->bounds);
    leafSAH /= area(bvh->bounds);
    assert(depth <= BVH4i::maxDepth);
  }

  std::string BVH4iStatistics::str()  
  {
    std::ostringstream stream;
    size_t bytesNodes = numNodes*sizeof(Node);
    size_t bytesTris  = numPrimBlocks*bvh->primTy.bytes;
    //size_t numVertices = bvh->numVertices;
    //size_t bytesVertices = numVertices*sizeof(Vec3fa); 
    size_t bytesTotal = bytesNodes+bytesTris;//+bytesVertices;
    size_t bytesTotalAllocated = bvh->bytes();
    double leafFill1 = double(numPrims)/double(bvh->primTy.blockSize*numPrimBlocks);
    double leafFill4 = double(numPrims)/double(4.0*numPrimBlocks4);
    stream.setf(std::ios::fixed, std::ios::floatfield);
    stream.precision(4);
    stream << "  sah = " << bvhSAH << ", leafSAH = " << leafSAH;
    stream.setf(std::ios::fixed, std::ios::floatfield);
    stream.precision(1);
    stream << ", depth = " << depth;
    stream << ", size = " << bytesTotalAllocated/1E6 << " MB (" << bytesTotal/1E6 << " MB used)" << std::endl;
    stream.precision(1);
    stream << "  nodes = "  << numNodes << " "
           << "(" << bytesNodes/1E6  << " MB) "
           << "(" << 100.0*double(bytesNodes)/double(bytesTotal) << "% of total) "
           << "(" << 100.0*(numNodes-1+numLeaves)/double(BVH4i::N*numNodes) << "% used)" 
           << std::endl;
    stream << "  leaves = " << numLeaves << " "
           << "(" << bytesTris/1E6  << " MB) "
           << "(" << 100.0*double(bytesTris)/double(bytesTotal) << "% of total) "
           << "(" << 100.0*leafFill1 << "% used, " << 100.0*leafFill4 << "% used)" 
           << std::endl;
    //stream << "  vertices = " << numVertices << " "
    //       << "(" << bytesVertices/1E6 << " MB) " 
    //       << "(" << 100.0*double(bytesVertices)/double(bytesTotal) << "% of total) "
    //       << "(" << 100.0*12.0f/float(sizeof(Vec3fa)) << "% used)" 
    //       << std::endl;
    return stream.str();
  }

  void BVH4iStatistics::statistics(NodeRef node, BBox3f bounds, size_t& depth)
  {
    float A = bounds.empty() ? 0.0f : area(bounds);
    
    if (node.isNode())
    {
      numNodes++;
      depth = 0;
      size_t cdepth = 0;
      Node* n = node.node(bvh->nodePtr());

      bvhSAH += A*BVH4i::travCost;
      for (size_t i=0; i<BVH4i::N; i++) {
        statistics(n->child(i),n->bounds(i),cdepth); 
        depth=max(depth,cdepth);
      }
      for (size_t i=0; i<BVH4i::N; i++) {
        if (n->child(i) == BVH4i::emptyNode) {
          for (; i<BVH4i::N; i++) {
            if (n->child(i) != BVH4i::emptyNode)
              throw std::runtime_error("invalid node");
          }
          break;
        }
      }    
      depth++;
      return;
    }
    else
    {
      depth = 0;
      unsigned int num; const char* tri = node.leaf(bvh->triPtr(),num);
      if (!num) return;
      
      numLeaves++;
      numPrimBlocks += num;
      size_t prims = 0;
      for (size_t i=0; i<num; i++) {
       prims += bvh->primTy.size(tri+i*bvh->primTy.bytes);
      }
      numPrims += prims;
      numPrimBlocks4 += (prims+3)/4;
      float sah = A * bvh->primTy.intCost * num;
      bvhSAH += sah;
      leafSAH += sah;
    }
  }
}

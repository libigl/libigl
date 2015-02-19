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

#include "bvh8_statistics.h"

namespace embree
{
  BVH8Statistics::BVH8Statistics (BVH8* bvh) : bvh(bvh)
  {
    numNodes = numLeaves = numPrimBlocks = numPrims = depth = 0;
    bvhSAH = leafSAH = 0.0f;
    statistics(bvh->root,bvh->bounds,depth);
    bvhSAH /= area(bvh->bounds);
    leafSAH /= area(bvh->bounds);
    assert(depth <= BVH8::maxDepth);
  }

  size_t BVH8Statistics::bytesUsed()
  {
    size_t bytesNodes = numNodes*sizeof(Node);
    size_t bytesTris  = numPrimBlocks*bvh->primTy.bytes;
    size_t numVertices = bvh->numVertices;
    size_t bytesVertices = numVertices*sizeof(Vec3fa); 
    return bytesNodes+bytesTris+bytesVertices;
  }

  std::string BVH8Statistics::str()  
  {
    std::ostringstream stream;
    size_t bytesNodes = numNodes*sizeof(Node);
    size_t bytesTris  = numPrimBlocks*bvh->primTy.bytes;
    size_t numVertices = bvh->numVertices;
    size_t bytesVertices = numVertices*sizeof(Vec3fa); 
    size_t bytesTotal = bytesNodes+bytesTris+bytesVertices;
    size_t bytesTotalAllocated = bvh->bytesAllocated();
    stream.setf(std::ios::fixed, std::ios::floatfield);
    stream << "  primitives = " << bvh->numPrimitives << ", vertices = " << bvh->numVertices << std::endl;
    stream.setf(std::ios::scientific, std::ios::floatfield);
    stream.precision(4);
    stream << "  sah = " << bvhSAH << ", leafSAH = " << leafSAH;
    stream.setf(std::ios::fixed, std::ios::floatfield);
    stream.precision(1);
    stream << ", depth = " << depth << std::endl;
    stream << "  used = " << bytesTotal/1E6 << " MB, allocated = " << bytesTotalAllocated/1E6 << " MB, perPrimitive = " << double(bytesTotal)/double(bvh->numPrimitives) << " B" << std::endl;
    stream.precision(1);
    stream << "  nodes = "  << numNodes << " "
           << "(" << bytesNodes/1E6  << " MB) "
           << "(" << 100.0*double(bytesNodes)/double(bytesTotal) << "% of total) "
           << "(" << 100.0*(numNodes-1+numLeaves)/(BVH8::N*numNodes) << "% used)" 
           << std::endl;
    stream << "  leaves = " << numLeaves << " "
           << "(" << bytesTris/1E6  << " MB) "
           << "(" << 100.0*double(bytesTris)/double(bytesTotal) << "% of total) "
           << "(" << 100.0*double(numPrims)/double(bvh->primTy.blockSize*numPrimBlocks) << "% used)" 
           << std::endl;
    stream << "  vertices = " << numVertices << " "
           << "(" << bytesVertices/1E6 << " MB) " 
           << "(" << 100.0*double(bytesVertices)/double(bytesTotal) << "% of total) "
           << "(" << 100.0*12.0f/float(sizeof(Vec3fa)) << "% used)" 
           << std::endl;
    return stream.str();
  }

  void BVH8Statistics::statistics(NodeRef node, const BBox3fa& bounds, size_t& depth)
  {
    float A = bounds.empty() ? 0.0f : area(bounds);

    if (node.isNode())
    {
      numNodes++;
      depth = 0;
      size_t cdepth = 0;
      Node* n = node.node();
      bvhSAH += A*BVH8::travCost;
      for (size_t i=0; i<BVH8::N; i++) {
        statistics(n->child(i),n->bounds(i),cdepth); 
        depth=max(depth,cdepth);
      }
      for (size_t i=0; i<BVH8::N; i++) {
        if (n->child(i) == BVH8::emptyNode) {
          for (; i<BVH8::N; i++) {
            if (n->child(i) != BVH8::emptyNode)
	      THROW_RUNTIME_ERROR("invalid node");
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
      size_t num; const char* tri = node.leaf(num);
      if (!num) return;
      
      numLeaves++;
      numPrimBlocks += num;
      for (size_t i=0; i<num; i++) {
        numPrims += bvh->primTy.size(tri+i*bvh->primTy.bytes);
      }
      float sah = A * bvh->primTy.intCost * num;
      bvhSAH += sah;
      leafSAH += sah;
    }
  }
}

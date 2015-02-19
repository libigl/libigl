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

#include "bvh4_statistics.h"

namespace embree
{
  BVH4Statistics::BVH4Statistics (BVH4* bvh) : bvh(bvh)
  {
    numAlignedNodes = numUnalignedNodes = 0;
    numAlignedNodesMB = numUnalignedNodesMB = 0;
    numLeaves = numPrims = depth = 0;
    childrenAlignedNodes = childrenUnalignedNodes = 0;
    childrenAlignedNodesMB = childrenUnalignedNodesMB = 0;
    bvhSAH = 0.0f;
    hash = 0;
    float A = max(0.0f,halfArea(bvh->bounds));
    statistics(bvh->root,A,depth);
    bvhSAH /= area(bvh->bounds);
    assert(depth <= BVH4::maxDepth);
  }

  size_t BVH4Statistics::bytesUsed() const
  {
    size_t bytesAlignedNodes = numAlignedNodes*sizeof(AlignedNode);
    size_t bytesUnalignedNodes = numUnalignedNodes*sizeof(UnalignedNode);
    size_t bytesAlignedNodesMB = numAlignedNodesMB*sizeof(BVH4::NodeMB);
    size_t bytesUnalignedNodesMB = numUnalignedNodesMB*sizeof(BVH4::UnalignedNodeMB);
    size_t bytesPrims  = numPrims*bvh->primTy.bytes;
    size_t numVertices = bvh->numVertices;
    size_t bytesVertices = numVertices*sizeof(Vec3fa); 
    return bytesAlignedNodes+bytesUnalignedNodes+bytesAlignedNodesMB+bytesUnalignedNodesMB+bytesPrims+bytesVertices;
  }

  std::string BVH4Statistics::str()  
  {
    std::ostringstream stream;
    size_t bytesAlignedNodes = numAlignedNodes*sizeof(AlignedNode);
    size_t bytesUnalignedNodes = numUnalignedNodes*sizeof(UnalignedNode);
    size_t bytesAlignedNodesMB = numAlignedNodesMB*sizeof(BVH4::NodeMB);
    size_t bytesUnalignedNodesMB = numUnalignedNodesMB*sizeof(BVH4::UnalignedNodeMB);
    size_t bytesPrims  = numPrims*bvh->primTy.bytes;
    size_t numVertices = bvh->numVertices;
    size_t bytesVertices = numVertices*sizeof(Vec3fa); 
    size_t bytesTotal = bytesAlignedNodes+bytesUnalignedNodes+bytesAlignedNodesMB+bytesUnalignedNodesMB+bytesPrims+bytesVertices;
    //size_t bytesTotalAllocated = bvh->alloc.bytes();
    stream.setf(std::ios::fixed, std::ios::floatfield);
    stream << "  primitives = " << bvh->numPrimitives << ", vertices = " << bvh->numVertices << ", hash= " << hash << std::endl;
    stream.precision(4);
    stream << "  sah = " << bvhSAH;
    stream.setf(std::ios::fixed, std::ios::floatfield);
    stream.precision(1);
    stream << ", depth = " << depth << std::endl;
    stream << "  used = " << bytesTotal/1E6 << " MB, perPrimitive = " << double(bytesTotal)/double(bvh->numPrimitives) << " B" << std::endl;
    stream.precision(1);
    if (numAlignedNodes) {
      stream << "  alignedNodes = "  << numAlignedNodes << " "
	     << "(" << 100.0*double(childrenAlignedNodes)/double(BVH4::N*numAlignedNodes) << "% filled) " 
	     << "(" << bytesAlignedNodes/1E6  << " MB) " 
	     << "(" << 100.0*double(bytesAlignedNodes)/double(bytesTotal) << "% of total)"
	     << std::endl;
    }
    if (numUnalignedNodes) {
      stream << "  unalignedNodes = "  << numUnalignedNodes << " "
	     << "(" << 100.0*double(childrenUnalignedNodes)/double(BVH4::N*numUnalignedNodes) << "% filled) " 
	     << "(" << bytesUnalignedNodes/1E6  << " MB) " 
	     << "(" << 100.0*double(bytesUnalignedNodes)/double(bytesTotal) << "% of total)"
	     << std::endl;
    }
    if (numAlignedNodesMB) {
      stream << "  alignedNodesMB = "  << numAlignedNodesMB << " "
	     << "(" << 100.0*double(childrenAlignedNodesMB)/double(BVH4::N*numAlignedNodesMB) << "% filled) " 
	     << "(" << bytesAlignedNodesMB/1E6  << " MB) " 
	     << "(" << 100.0*double(bytesAlignedNodesMB)/double(bytesTotal) << "% of total)"
	     << std::endl;
    }
    if (numUnalignedNodesMB) {
      stream << "  unalignedNodesMB = "  << numUnalignedNodesMB << " "
	     << "(" << 100.0*double(childrenUnalignedNodesMB)/double(BVH4::N*numUnalignedNodesMB) << "% filled) " 
	     << "(" << bytesUnalignedNodesMB/1E6  << " MB) " 
	     << "(" << 100.0*double(bytesUnalignedNodesMB)/double(bytesTotal) << "% of total)"
	     << std::endl;
    }
    stream << "  leaves = " << numLeaves << " "
           << "(" << bytesPrims/1E6  << " MB) "
           << "(" << 100.0*double(bytesPrims)/double(bytesTotal) << "% of total)"
           << std::endl;
    stream << "  vertices = " << numVertices << " "
           << "(" << bytesVertices/1E6 << " MB) " 
           << "(" << 100.0*double(bytesVertices)/double(bytesTotal) << "% of total) "
           << "(" << 100.0*12.0f/float(sizeof(Vec3fa)) << "% used)" 
           << std::endl;
    return stream.str();
  }

  void BVH4Statistics::statistics(NodeRef node, const float A, size_t& depth)
  {
    if (node.isNode())
    {
      hash += 0x1234;
      numAlignedNodes++;
      AlignedNode* n = node.node();
      bvhSAH += A*BVH4::travCostAligned;

      depth = 0;
      for (size_t i=0; i<BVH4::N; i++) {
        if (n->child(i) != BVH4::emptyNode) childrenAlignedNodes++;
        const float Ai = max(0.0f,halfArea(n->extend(i)));
        size_t cdepth; statistics(n->child(i),Ai,cdepth); 
        depth=max(depth,cdepth);
      }
      depth++;
      hash += 0x76767*depth;
    }
    else if (node.isUnalignedNode())
    {
      hash += 0x1232344;
      numUnalignedNodes++;
      UnalignedNode* n = node.unalignedNode();
      bvhSAH += A*BVH4::travCostUnaligned;

      depth = 0;
      for (size_t i=0; i<BVH4::N; i++) {
        if (n->child(i) != BVH4::emptyNode) childrenUnalignedNodes++;
        const float Ai = max(0.0f,halfArea(n->extend(i)));
        size_t cdepth; statistics(n->child(i),Ai,cdepth); 
        depth=max(depth,cdepth);
      }
      depth++;
      hash += 0x76767*depth;
    }
    else if (node.isNodeMB())
    {
      hash += 0xEF343;
      numAlignedNodesMB++;
      BVH4::NodeMB* n = node.nodeMB();
      bvhSAH += A*BVH4::travCostAligned;

      depth = 0;
      for (size_t i=0; i<BVH4::N; i++) {
        if (n->child(i) != BVH4::emptyNode) childrenAlignedNodesMB++;
        const float Ai = max(0.0f,halfArea(n->extend0(i)));
        size_t cdepth; statistics(n->child(i),Ai,cdepth); 
        depth=max(depth,cdepth);
      }
      depth++;
      hash += 0x76767*depth;
    }
    else if (node.isUnalignedNodeMB())
    {
      hash += 0x1EEF4;
      numUnalignedNodesMB++;
      BVH4::UnalignedNodeMB* n = node.unalignedNodeMB();
      bvhSAH += A*BVH4::travCostUnaligned;

      depth = 0;
      for (size_t i=0; i<BVH4::N; i++) {
        if (n->child(i) != BVH4::emptyNode) childrenUnalignedNodesMB++;
        const float Ai = max(0.0f,halfArea(n->extend0(i)));
        size_t cdepth; statistics(n->child(i),Ai,cdepth); 
        depth=max(depth,cdepth);
      }
      depth++;
      hash += 0x76767*depth;
    }
    else
    {
      depth = 0;
      size_t num; const char* tri = node.leaf(num);
      hash += 0xDD776*num+0x878;
      if (!num) return;

      hash += bvh->primTy.hash(tri,num);
      
      numLeaves++;
      numPrims += num;
      float sah = A * BVH4::intCost * num;
      bvhSAH += sah;
    }
  }
}

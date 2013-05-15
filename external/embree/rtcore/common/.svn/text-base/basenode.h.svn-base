// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#ifndef __EMBREE_BASE_NODE_H__
#define __EMBREE_BASE_NODE_H__

#include "../common/accel.h"

namespace embree
{
  /*! Pointer that points to a Node or a list of Triangles */
  template<typename Node, int alignment>
    struct BaseNode
  {
    /*! Masks the bits that store the number of items per leaf. */
    static const size_t mask =  (1 << (alignment-1))-1;  

    /*! Maximal number of triangle blocks per leaf. */
    static const size_t maxLeafBlocks = mask-1;
    
    /*! Empty node */
    static const size_t empty = 1;

    /*! checks if this is an empty leaf */
    __forceinline int isEmptyLeaf() const { return ((size_t)this & (size_t)mask) == empty; }

    /*! Checks if this is an barrier. A barrier tells the top level tree rotation of how deep to enter the tree. */
    __forceinline int isBarrier() const { return ((size_t)this >> (size_t)(alignment-1)) & 1; }

    /*! Clears the barrier bit. */
    __forceinline BaseNode* setBarrier() const { return (BaseNode*)((size_t)this | (size_t)(1 << (alignment-1))); }

    /*! Clears the barrier bit. */
    __forceinline BaseNode* clearBarrier() const { return (BaseNode*)((size_t)this & ~(size_t)(1 << (alignment-1))); }

    /*! checks if this is a leaf */
    __forceinline int isLeaf() const { return ((size_t)this & (size_t)mask) != 0; }
    
    /*! checks if this is a node */
    __forceinline int isNode() const { return ((size_t)this & (size_t)mask) == 0; }
    
    /*! returns node pointer */
    __forceinline Node* node() const { 
      assert(isNode());
      return (Node*)this;
    }
    
    /*! returns leaf pointer */
    __forceinline char* leaf(size_t& num) const {
      assert(isLeaf());
      num = ((size_t)this & (size_t)mask)-1;
      return (char*)((size_t)this & ~(size_t)mask);
    }
    
    /*! encodes a node */
    __forceinline static BaseNode* encodeNode(Node* node) { 
      return (BaseNode*)node;
    }
    
    /*! encodes a leaf */
    __forceinline static BaseNode* encodeLeaf(char* tri, size_t num) {
      assert(!((size_t)tri & (size_t)mask));
#if defined(_DEBUG)
      if (num > maxLeafBlocks) throw std::runtime_error("ERROR: Loosing triangles during build.");
#else
      if (num > maxLeafBlocks) std::cerr << "WARNING: Loosing triangles during build." << std::endl;
#endif
      return (BaseNode*)((size_t)tri | (1+min(num,size_t(maxLeafBlocks))));
    }
  };

  /*! Computes half surface area of box. */
  __forceinline float halfArea3f(const BBox<ssef>& box) {
    const ssef d = size(box);
    const ssef a = d*shuffle<1,2,0,3>(d);
    return a[0]+a[1]+a[2];
  }
}

#endif

  

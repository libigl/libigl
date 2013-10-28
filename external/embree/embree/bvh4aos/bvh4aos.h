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

#ifndef __EMBREE_BVH4AOS_H__
#define __EMBREE_BVH4AOS_H__

#include "common/accel.h"
#include "common/alloc.h"
#include "geometry/triangle.h"
#include "simd/mic.h"

namespace embree
{
  /*! Multi BVH with 4 children. Each node stores the bounding box of
   * it's 4 children in AOS format as well as a 4 child indices. */
  class BVH4AOS : public Accel
  {
  public:

    /*! forward declaration of node type */
    struct Node;

    /*! branching width of the tree */
    static const size_t N = 4;

    /*! Maximal depth of the BVH. */
    static const size_t maxBuildDepth = 32;
    static const size_t maxBuildDepthLeaves = maxBuildDepth+16;
    static const size_t maxDepth = maxBuildDepth+maxBuildDepthLeaves; // makes tree rotations of top part of tree safe
         
    /*! Each reference to a node consists of offset (26 bit), 
        is-leaf bit (1 bit), and number of triangle count (5 bits). */
    static const unsigned index_shift = 6;           
    static const unsigned leaf_mask = (1 << 5);
    static const unsigned items_mask = (1 << 5)-1;
    static const unsigned offset_mask = ~((1 << 6) - 1);
    static const size_t maxLeafBlocks = items_mask; 
    static const unsigned invalid_node = leaf_mask | (1 << 4);
    static const mic_f initial_node;

#define BVH4AOS_INVALID_NODE ((1 << 5) | (1 << 4))
#define BVH4AOS_LEAF_MASK ((unsigned)1 << 5)
#define BVH4AOS_OFFSET_MASK (~((1 << 6) - 1))
#define BVH4AOS_ITEMS_MASK ((1 << 5)-1)

    /*! Cost of one traversal step. */
    static const int travCost = 1;   

    static __forceinline unsigned isLeaf(unsigned id) { return (id & ((1 << 5))); }

    static __forceinline unsigned itemOffset(const unsigned children) {
      return children & BVH4AOS_OFFSET_MASK; // 6 bits instead of 7
    }
    
    static __forceinline unsigned itemCount(const unsigned children) {
      return children & BVH4AOS_ITEMS_MASK; // 6 bits instead of 7
    }

    /*! Reference to a node */
    struct NodeRef
    {
      /*! Default constructor */
      __forceinline NodeRef () {}

      /*! Construction from integer */
      __forceinline NodeRef (unsigned id) : id(id) { }

      /*! Cast to int */
      __forceinline operator unsigned() const { return id; }

      /*! set barrier bit */
      __forceinline void setBarrier() { this->id |= (1<<4); }

      /*! clear barrier bit */
      __forceinline void clearBarrier() { this->id &= ~(1<<4); }
      
      /*! test barrier bit */
      __forceinline bool isBarrier() const { return this->id & (1<<4); }

      /*! tests if node is invalid */
      __forceinline int isInvalid() const { return id == BVH4AOS::invalid_node; }

      /*! checks if this is a leaf */
      __forceinline int isLeaf() const { return (id & ((1 << 5))); }
      //__forceinline int isLeaf() const { return id < 0; }
      
      /*! checks if this is a node */
      __forceinline int isNode() const { return ((~id) & leaf_mask); }
      //__forceinline int isNode() const { return id >= 0; }
      
      /*! returns node pointer */
      __forceinline       Node* node(      void* base) const { assert(isNode()); return (      Node*)((      char*)base + id); }
      __forceinline const Node* node(const void* base) const { assert(isNode()); return (const Node*)((const char*)base + id); }
      
      /*! returns leaf pointer */
      __forceinline const char* leaf(const void* base, size_t& num) const {
        assert(isLeaf());
        num = id & items_mask;
        return (const char*)base + (id & offset_mask);
      }
      
      /*! encodes a node */
      __forceinline static NodeRef encodeNode(void* base, Node* node) { 
        ssize_t ofs = (size_t)node-(size_t)base;
        assert((ofs & ~offset_mask) == 0);
        return NodeRef((unsigned)ofs);
      }
      
      /*! encodes a leaf */
      __forceinline static NodeRef encodeLeaf(void* base, char* tri, size_t num) {
        ssize_t ofs = (size_t)tri-(size_t)base;
        assert((ofs & ~offset_mask) == 0);
#if defined(_DEBUG)
        if (num > maxLeafBlocks) throw std::runtime_error("ERROR: Loosing triangles during build.");
#else
        if (num > maxLeafBlocks) std::cerr << "WARNING: Loosing triangles during build." << std::endl;
#endif
        return NodeRef(leaf_mask | ofs | min(num,size_t(maxLeafBlocks)));
      }
      
    public:
      unsigned id;
    };

    /*! Node structure */
    struct Node
    {
      /*! Clears the node. */
      __forceinline void clear() { 
        (mic_f&) lower = initial_node;
        (mic_f&) upper = initial_node;
      }

      /*! Sets bounding box of child. */
      __forceinline void set(size_t i, const BBox3f& bounds) {
        lower[i].x = bounds.lower.x; lower[i].y = bounds.lower.y; lower[i].z = bounds.lower.z;
        upper[i].x = bounds.upper.x; upper[i].y = bounds.upper.y; upper[i].z = bounds.upper.z;
      }

      /*! Sets bounding box and ID of child. */
      __forceinline void set(size_t i, const BBox3f& bounds, NodeRef childID) { 
        set(i,bounds);
        lower[i].child = childID;
      }

      /*! Returns bounds of node. */
      __forceinline BBox3f bounds() const {
        return merge(bounds(0),bounds(1),bounds(2),bounds(3));
      }

      /*! Returns bounds of specified child. */
      __forceinline BBox3f bounds(size_t i) const 
      {
        BBox3f b(Vector3f(load4f(&lower[i])),Vector3f(load4f(&upper[i])));
        if (unlikely(b.lower == b.upper)) return BBox3f(empty);
        return b;
      }

      /*! Returns reference to specified child */
      __forceinline       NodeRef& child(size_t i)       { return lower[i].child; }
      __forceinline const NodeRef& child(size_t i) const { return lower[i].child; }

    public:
      struct NodeStruct {
        float x,y,z;           // x,y, and z coordinates of bounds
        NodeRef child;         // encodes 1 is-leaf bit, 25 offset bits, and 6 num-items bits
      } lower[4], upper[4];    // lower and upper bounds of all 4 children
    };

    /*! swap the children of two nodes */
    __forceinline static void swap(Node* a, size_t i, Node* b, size_t j) {
      std::swap(a->lower[i],b->lower[j]);
      std::swap(a->upper[i],b->upper[j]);
    }

  public:

    /*! BVH4 default constructor. */
    BVH4AOS (RTCGeometry* geom, 
	     const TriangleType& trity, 
             size_t bytesNodeArray     = 1LL << 32, 
	     size_t bytesTriangleArray = 1LL << 32)
      : Accel(geom,trity), maxLeafTris(maxLeafBlocks*trity.blockSize), root(NULL)
    {
      alloc_nodes = new LinearAllocatorPerThread;
      alloc_tris  = new LinearAllocatorPerThread;
    }

    /*! create function for registry */
    static Accel* create (RTCGeometry* geom, const TriangleType& trity) {
      return new BVH4AOS(geom,trity);
    }

    /*! initializes the acceleration structure */
    void init (size_t numNodes, size_t numTriangles);

    /*! clears the acceleration structure */
    void resize (size_t bytesNodeArray, size_t bytesTriangleArray);

    /*! return name of this acceleration structure */
    const std::string name() const { return "bvh4aos."+trity.name; }

    /*! prints statistics */
    void print();

    void clearBarrier(NodeRef& node)
    {
      if (node.isBarrier()) {
        node.clearBarrier();
        return;
      }
      else if (node.isLeaf()) 
        return;
      else {
        Node* n = node.node(nodePtr());
        for (size_t c=0; c<4; c++)
          clearBarrier(n->child(c));
      }
    }

    /*! Data of the BVH */
  public:
    const size_t maxLeafTris;          //!< maximal number of triangles per leaf
    NodeRef root;                      //!< Root node (can also be a leaf).
    
    /*! Memory allocation */
  public:
    Ref<LinearAllocatorPerThread> alloc_nodes;
    Ref<LinearAllocatorPerThread> alloc_tris;

    __forceinline Node* allocNode(size_t thread) {
      Node* node = (Node*) alloc_nodes->malloc(thread,sizeof(Node),1 << index_shift); node->clear(); return node;
    }

    __forceinline char* allocTris(size_t thread, size_t num) {
      return (char*) alloc_tris->malloc(thread,num*trity.bytes,1 << index_shift);
    }

    __forceinline       void* nodePtr()       { return alloc_nodes->base(); }
    __forceinline const void* nodePtr() const { return alloc_nodes->base(); }

    __forceinline       void* triPtr()       { return alloc_tris->base(); }
    __forceinline const void* triPtr() const { return alloc_tris->base(); }
  };
}

#endif

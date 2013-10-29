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

#ifndef __EMBREE_BVH4I_H__
#define __EMBREE_BVH4I_H__

#include "../common/accel.h"
#include "../common/alloc.h"
#include "../geometry/triangle.h"

namespace embree
{
  /*! Multi BVH with 4 children. Each node stores the bounding box of
   * it's 4 children as well as a 4 child indices. */
  class BVH4i : public Accel
  {
  public:

    /*! forward declaration of node type */
    struct Node;

    /*! branching width of the tree */
    static const size_t N = 4;

    /*! Number of address bits the Node and Triangle are aligned
        to. Maximally 2^alignment-2 many triangle blocks per leaf are
        supported. */
    static const size_t alignment = 4;

    /*! Masks the bits that store the number of items per leaf. */
    static const unsigned offset_mask = 0x7FFFFFFF & (0xFFFFFFFF << alignment);
    static const unsigned barrier_mask = 1 << (alignment-1);  
    static const unsigned items_mask = barrier_mask-1;  
    
    /*! Empty node */
    static const unsigned emptyNode = 0x80000000;
      
    /*! Maximal depth of the BVH. */
    static const size_t maxBuildDepth = 32;
    static const size_t maxBuildDepthLeaves = maxBuildDepth+16;
    static const size_t maxDepth = maxBuildDepth+maxBuildDepthLeaves; // makes tree rotations of top part of tree safe
    
    /*! Maximal number of triangle blocks in a leaf. */
    static const size_t maxLeafBlocks = items_mask-1;

    /*! Cost of one traversal step. */
    static const int travCost = 1;

    /*! References a Node or list of Triangles */
    struct NodeRef
    {
      /*! Default constructor */
      __forceinline NodeRef () {}

      /*! Construction from integer */
      __forceinline NodeRef (unsigned id) : id(id) { }

      /*! Cast to unsigned */
      __forceinline operator unsigned() const { return id; }

      /*! Clears the barrier bit. */
      __forceinline void setBarrier() { id |= barrier_mask; }
      
      /*! Clears the barrier bit. */
      __forceinline void clearBarrier() { id &= ~barrier_mask; }

      /*! Checks if this is an barrier. A barrier tells the top level tree rotation of how deep to enter the tree. */
      __forceinline unsigned isBarrier() const { return id & barrier_mask; }
     
      /*! checks if this is a leaf */
      __forceinline unsigned isLeaf() const { return (int) id < 0; }
      
      /*! checks if this is a node */
      __forceinline unsigned isNode() const { return (int) id >= 0; }
      
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
      __forceinline static NodeRef encodeNode(void* base, Node* node) 
      { 
        ssize_t ofs = (size_t)node-(size_t)base;
        assert(ofs >= 0);
        assert((ofs & ~(size_t)offset_mask) == 0);
        return NodeRef(ofs);
      }
      
      /*! encodes a leaf */
      __forceinline static NodeRef encodeLeaf(void* base, char* tri, size_t num) 
      {
        ssize_t ofs = (size_t)tri-(size_t)base;
        assert(ofs >= 0);
        assert((ofs & ~(size_t)offset_mask) == 0);
#if defined(_DEBUG)
        if (num > (size_t)maxLeafBlocks) throw std::runtime_error("ERROR: Loosing triangles during build.");
#else
        if (num > (size_t)maxLeafBlocks) std::cout << "WARNING: Loosing triangles during build." << std::endl;
#endif
        return NodeRef(0x80000000 | ofs | min(num,(size_t)maxLeafBlocks));
      }

    private:
      unsigned id;
    };

    /*! BVH4 Node */
    struct Node
    {
      /*! Clears the node. */
      __forceinline void clear()  {
        lower_x = lower_y = lower_z = ssef(1E10); 
        upper_x = upper_y = upper_z = ssef(1E10);
        children[0] = children[1] = children[2] = children[3] = emptyNode;
      }

      /*! Sets bounding box of child. */
      __forceinline void set(size_t i, const BBox3f& bounds) {
        lower_x[i] = bounds.lower.x; lower_y[i] = bounds.lower.y; lower_z[i] = bounds.lower.z;
        upper_x[i] = bounds.upper.x; upper_y[i] = bounds.upper.y; upper_z[i] = bounds.upper.z;
      }

      /*! Sets bounding box and ID of child. */
      __forceinline void set(size_t i, const BBox3f& bounds, const NodeRef& childID) {
        set(i,bounds);
        children[i] = childID;
      }

      /*! Returns bounds of node. */
      __forceinline BBox3f bounds() const {
        sseb invalid = (lower_x == upper_x) & (lower_y == upper_y) & (lower_z == upper_z);
        ssef lx = select(invalid,pos_inf,lower_x);
        ssef ly = select(invalid,pos_inf,lower_y);
        ssef lz = select(invalid,pos_inf,lower_z);
        ssef ux = select(invalid,neg_inf,upper_x);
        ssef uy = select(invalid,neg_inf,upper_y);
        ssef uz = select(invalid,neg_inf,upper_z);
        return BBox3f(Vector3f(reduce_min(lx),reduce_min(ly),reduce_min(lz)),
                      Vector3f(reduce_max(ux),reduce_max(uy),reduce_max(uz)));
      }

      /*! Returns bounds of specified child. */
      __forceinline BBox3f bounds(size_t i) const {
        BBox3f box(Vector3f(lower_x[i],lower_y[i],lower_z[i]),Vector3f(upper_x[i],upper_y[i],upper_z[i]));
        if (likely(box.lower != box.upper)) return box;
        return BBox3f(empty);
      }

      /*! Returns reference to specified child */
      __forceinline       NodeRef& child(size_t i)       { return children[i]; }
      __forceinline const NodeRef& child(size_t i) const { return children[i]; }

    public:
      ssef_m lower_x;           //!< X dimension of lower bounds of all 4 children.
      ssef_m upper_x;           //!< X dimension of upper bounds of all 4 children.
      ssef_m lower_y;           //!< Y dimension of lower bounds of all 4 children.
      ssef_m upper_y;           //!< Y dimension of upper bounds of all 4 children.
      ssef_m lower_z;           //!< Z dimension of lower bounds of all 4 children.
      ssef_m upper_z;           //!< Z dimension of upper bounds of all 4 children.
      NodeRef children[4];      //!< Pointer to the 4 children (can be a node or leaf)
    };

    /*! swap the children of two nodes */
    __forceinline static void swap(Node* a, size_t i, Node* b, size_t j)
    {
      std::swap(a->children[i],b->children[j]);
      std::swap(a->lower_x[i],b->lower_x[j]);
      std::swap(a->lower_y[i],b->lower_y[j]);
      std::swap(a->lower_z[i],b->lower_z[j]);
      std::swap(a->upper_x[i],b->upper_x[j]);
      std::swap(a->upper_y[i],b->upper_y[j]);
      std::swap(a->upper_z[i],b->upper_z[j]);
    }
    
  public:

    /*! BVH4 default constructor. */
    BVH4i (RTCGeometry* geom, const TriangleType& trity)
      : Accel(geom,trity), maxLeafTris(maxLeafBlocks*trity.blockSize), root(emptyNode)
    {
      alloc_nodes = new LinearAllocatorPerThread;
      alloc_tris  = new LinearAllocatorPerThread;
    }

    /*! create function for registry */
    static Accel* create (RTCGeometry* geom, const TriangleType& trity) {
      return new BVH4i(geom,trity);
    }

    /*! initializes the acceleration structure */
    void init (size_t numNodes, size_t numTriangles);

    /*! return name of this acceleration structure */
    const std::string name() const { return "bvh4i."+trity.name; }

    /*! prints statistics */
    void print();

    /*! Clears the barrier bits. */
    void clearBarrier(NodeRef& node)
    {
      if (node.isBarrier()) 
        node.clearBarrier();
      else if (!node.isLeaf()) {
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
      Node* node = (Node*) alloc_nodes->malloc(thread,sizeof(Node),1 << alignment); node->clear(); return node;
    }

    __forceinline char* allocTris(size_t thread, size_t num) {
      return (char*) alloc_tris->malloc(thread,num*trity.bytes,1 << alignment);
    }

    __forceinline       void* nodePtr()       { return alloc_nodes->base(); }
    __forceinline const void* nodePtr() const { return alloc_nodes->base(); }

    __forceinline       void* triPtr()       { return alloc_tris->base(); }
    __forceinline const void* triPtr() const { return alloc_tris->base(); }
  };
}

#endif

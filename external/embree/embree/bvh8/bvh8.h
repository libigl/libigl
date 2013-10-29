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

#ifndef __EMBREE_BVH8_H__
#define __EMBREE_BVH8_H__

#include "../common/accel.h"
#include "../common/alloc.h"
#include "../geometry/triangle.h"

namespace embree
{
  /*! Multi BVH with 8 children. Each node stores the bounding box of
   * it's 8 children as well as a 8 child pointers. */
  class BVH8 : public Accel
  {
  public:
    
    /*! forward declaration of node type */
    struct Node;

    /*! branching width of the tree */
    static const size_t N = 8;

    /*! Number of address bits the Node and Triangle are aligned
        to. Maximally 2^alignment-2 many triangle blocks per leaf are
        supported. */
    static const size_t alignment = 5;

    /*! Masks the bits that store the number of items per leaf. */
    static const size_t align_mask = (1 << alignment)-1;  
    static const size_t items_mask = (1 << (alignment-1))-1;  

    /*! Empty node */
    static const size_t emptyNode = 1;
      
    /*! Maximal depth of the BVH. */
    static const size_t maxBuildDepth = 32;
    static const size_t maxBuildDepthLeaves = maxBuildDepth+16;
    static const size_t maxDepth = maxBuildDepth+maxBuildDepthLeaves; // makes tree rotations of top part of tree safe
    
    /*! Maximal number of triangle blocks in a leaf. */
    static const size_t maxLeafBlocks = items_mask-1;

    /*! Cost of one traversal step. */
    static const int travCost = 1;

    /*! Pointer that points to a Node or a list of Triangles */
    struct NodeRef
    {
      typedef size_t Ty;

      /*! Default constructor */
      __forceinline NodeRef () {}

      /*! Construction from integer */
      __forceinline NodeRef (size_t ptr) : ptr(ptr) { }

      /*! Cast to Ty */
      __forceinline operator Ty() const { return ptr; }

      /*! Clears the barrier bit. */
      __forceinline void setBarrier() { ptr |= (size_t)(1 << (alignment-1)); }
      
      /*! Clears the barrier bit. */
      __forceinline void clearBarrier() { ptr &= ~(size_t)(1 << (alignment-1)); }

      /*! Checks if this is an barrier. A barrier tells the top level tree rotation of how deep to enter the tree. */
      __forceinline int isBarrier() const { return (ptr >> (size_t)(alignment-1)) & 1; }

      /*! checks if this is a leaf */
      __forceinline int isLeaf() const { return (ptr & (size_t)align_mask) != 0; }
      
      /*! checks if this is a node */
      __forceinline int isNode() const { return (ptr & (size_t)align_mask) == 0; }
      
      /*! returns node pointer */
      __forceinline       Node* node(      void* base) const { assert(isNode()); return (      Node*)ptr; }
      __forceinline const Node* node(const void* base) const { assert(isNode()); return (const Node*)ptr; }
      
      /*! returns leaf pointer */
      __forceinline const char* leaf(const void* base, size_t& num) const {
        assert(isLeaf());
        num = (ptr & (size_t)items_mask)-1;
        return (char*)(ptr & ~(size_t)align_mask);
      }
      
      /*! encodes a node */
      __forceinline static NodeRef encodeNode(void* base, Node* node) { 
        return NodeRef((size_t) node);
      }
      
      /*! encodes a leaf */
      __forceinline static NodeRef encodeLeaf(void* base, char* tri, size_t num) {
        assert(!((size_t)tri & align_mask)); 
#if defined(_DEBUG)
        if (num > (size_t)maxLeafBlocks) throw std::runtime_error("ERROR: Loosing triangles during build.");
#else
        if (num > (size_t)maxLeafBlocks) std::cout << "WARNING: Loosing triangles during build." << std::endl;
#endif
        return NodeRef((size_t)tri | (1+min(num,(size_t)maxLeafBlocks)));
      }

    private:
      size_t ptr;
    };

    /*! BVH8 Node */
    struct Node
    {
      /*! Clears the node. */
      __forceinline void clear()  {
        lower_x = lower_y = lower_z = pos_inf; 
        upper_x = upper_y = upper_z = neg_inf;
        for (size_t i=0; i<N; i++) children[i] = emptyNode;
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
        return BBox3f(Vector3f(reduce_min(lower_x),reduce_min(lower_y),reduce_min(lower_z)),
                      Vector3f(reduce_max(upper_x),reduce_max(upper_y),reduce_max(upper_z)));
      }

      /*! Returns bounds of specified child. */
      __forceinline BBox3f bounds(size_t i) const {
        return BBox3f(Vector3f(lower_x[i],lower_y[i],lower_z[i]),
                      Vector3f(upper_x[i],upper_y[i],upper_z[i]));
      }

      /*! Returns reference to specified child */
      __forceinline       NodeRef& child(size_t i)       { return children[i]; }
      __forceinline const NodeRef& child(size_t i) const { return children[i]; }

    public:
      avxf lower_x;           //!< X dimension of lower bounds of all 8 children.
      avxf upper_x;           //!< X dimension of upper bounds of all 8 children.
      avxf lower_y;           //!< Y dimension of lower bounds of all 8 children.
      avxf upper_y;           //!< Y dimension of upper bounds of all 8 children.
      avxf lower_z;           //!< Z dimension of lower bounds of all 8 children.
      avxf upper_z;           //!< Z dimension of upper bounds of all 8 children.
      NodeRef children[N];      //!< Pointer to the 8 children (can be a node or leaf)
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

    /*! BVH8 default constructor. */
    BVH8 (RTCGeometry* geom, const TriangleType& trity)
      : Accel(geom,trity), maxLeafTris(maxLeafBlocks*trity.blockSize), root(emptyNode) {}
    
    /*! create function for registry */
    static Accel* create (RTCGeometry* geom, const TriangleType& trity) {
      return new BVH8(geom,trity);
    }

    /*! initializes the acceleration structure */
    void init (size_t numNodes, size_t numTriangles);

    /*! return name of this acceleration structure */
    const std::string name() const { return "bvh8."+trity.name; }

    /*! prints statistics */
    void print();

    /*! Clears the barrier bits. */
    void clearBarrier(NodeRef& node)
    {
      if (node.isBarrier())
        node.clearBarrier();
      else if (!node.isLeaf()) {
        Node* n = node.node(nodePtr());
        for (size_t c=0; c<8; c++)
          clearBarrier(n->child(c));
      }
    }

    /*! Data of the BVH */
  public:
    const size_t maxLeafTris;          //!< maximal number of triangles per leaf
    NodeRef root;                      //!< Root node (can also be a leaf).

  /*! Memory allocation */
  public:
    AllocatorPerThread alloc;          //!< allocator for nodes and triangles
 
    __forceinline Node* allocNode(size_t thread) {
      Node* node = (Node*) alloc.malloc(thread,sizeof(Node),1 << alignment); node->clear(); return node;
    }

    __forceinline char* allocTris(size_t thread, size_t num) {
      return (char*) alloc.malloc(thread,num*trity.bytes,1 << alignment);
    }
  };
}

#endif

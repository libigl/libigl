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

#ifndef __EMBREE_BVH4_H__
#define __EMBREE_BVH4_H__

#include "embree2/rtcore.h"
#include "common/alloc.h"
#include "common/accel.h"
#include "common/scene.h"
#include "geometry/primitive.h"

namespace embree
{
  /*! Multi BVH with 4 children. Each node stores the bounding box of
   * it's 4 children as well as 4 child pointers. */
  class BVH4 : public Bounded
  {
	  ALIGNED_CLASS;
  public:
    
    /*! forward declaration of node type */
    struct Node;

    /*! branching width of the tree */
    static const size_t N = 4;

    /*! Number of address bits the Node and primitives are aligned
        to. Maximally 2^alignment-2 many primitive blocks per leaf are
        supported. */
    static const size_t alignment = 4;

    /*! Masks the bits that store the number of items per leaf. */
    static const size_t align_mask = (1 << alignment)-1;  
    static const size_t items_mask = (1 << (alignment-1))-1;  

    /*! Empty node */
    static const size_t emptyNode = 1;

    /*! Invalid node, used as marker in traversal */
    static const size_t invalidNode = (((size_t)-1) & (~items_mask)) | 1;
      
    /*! Maximal depth of the BVH. */
    static const size_t maxBuildDepth = 32;
    static const size_t maxBuildDepthLeaf = maxBuildDepth+16;
    static const size_t maxDepth = maxBuildDepthLeaf+maxBuildDepthLeaf+maxBuildDepth;
    
    /*! Maximal number of primitive blocks in a leaf. */
    static const size_t maxLeafBlocks = items_mask-1;

    /*! Cost of one traversal step. */
    static const int travCost = 1;

    /*! Pointer that points to a node or a list of primitives */
    struct NodeRef
    {
      /*! Default constructor */
      __forceinline NodeRef () {}

      /*! Construction from integer */
      __forceinline NodeRef (size_t ptr) : ptr(ptr) { }

      /*! Cast to size_t */
      __forceinline operator size_t() const { return ptr; }

      /*! Sets the barrier bit. */
      __forceinline void setBarrier() { ptr |= (size_t)(1 << (alignment-1)); }
      
      /*! Clears the barrier bit. */
      __forceinline void clearBarrier() { ptr &= ~(size_t)(1 << (alignment-1)); }

      /*! Checks if this is an barrier. A barrier tells the top level tree rotations how deep to enter the tree. */
      __forceinline int isBarrier() const { return (ptr >> (size_t)(alignment-1)) & 1; }

      /*! checks if this is a leaf */
      __forceinline int isLeaf() const { return (ptr & (size_t)align_mask) != 0; }
      
      /*! checks if this is a node */
      __forceinline int isNode() const { return (ptr & (size_t)align_mask) == 0; }
      
      /*! returns node pointer */
      __forceinline       Node* node()       { assert(isNode()); return (      Node*)ptr; }
      __forceinline const Node* node() const { assert(isNode()); return (const Node*)ptr; }
      
      /*! returns leaf pointer */
      __forceinline char* leaf(size_t& num) const {
        assert(isLeaf());
        num = (ptr & (size_t)items_mask)-1;
        return (char*)(ptr & ~(size_t)align_mask);
      }

    private:
      size_t ptr;
    };

    /*! BVH4 Node */
    struct Node
    {
      /*! Clears the node. */
      __forceinline void clear() {
        lower_x = lower_y = lower_z = pos_inf; 
        upper_x = upper_y = upper_z = neg_inf;
        children[0] = children[1] = children[2] = children[3] = emptyNode;
      }

      /*! Sets bounding box of child. */
      __forceinline void set(size_t i, const BBox3f& bounds) 
      {
        assert(i < 4);
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
        const Vec3fa lower(reduce_min(lower_x),reduce_min(lower_y),reduce_min(lower_z));
        const Vec3fa upper(reduce_max(upper_x),reduce_max(upper_y),reduce_max(upper_z));
        return BBox3f(lower,upper);
      }

      /*! Returns bounds of specified child. */
      __forceinline BBox3f bounds(size_t i) const 
      {
        assert(i < 4);
        const Vec3fa lower(lower_x[i],lower_y[i],lower_z[i]);
        const Vec3fa upper(upper_x[i],upper_y[i],upper_z[i]);
        return BBox3f(lower,upper);
      }

      /*! Returns bounds of all children */
      __forceinline void bounds(BBox<ssef>& bounds0, BBox<ssef>& bounds1, BBox<ssef>& bounds2, BBox<ssef>& bounds3) const {
        transpose(lower_x,lower_y,lower_z,ssef(zero),bounds0.lower,bounds1.lower,bounds2.lower,bounds3.lower);
        transpose(upper_x,upper_y,upper_z,ssef(zero),bounds0.upper,bounds1.upper,bounds2.upper,bounds3.upper);
      }

      /*! Returns reference to specified child */
      __forceinline       NodeRef& child(size_t i)       { assert(i<4); return children[i]; }
      __forceinline const NodeRef& child(size_t i) const { assert(i<4); return children[i]; }

    public:
      ssef lower_x;           //!< X dimension of lower bounds of all 4 children.
      ssef upper_x;           //!< X dimension of upper bounds of all 4 children.
      ssef lower_y;           //!< Y dimension of lower bounds of all 4 children.
      ssef upper_y;           //!< Y dimension of upper bounds of all 4 children.
      ssef lower_z;           //!< Z dimension of lower bounds of all 4 children.
      ssef upper_z;           //!< Z dimension of upper bounds of all 4 children.
      NodeRef children[4];    //!< Pointer to the 4 children (can be a node or leaf)
    };

    /*! swap the children of two nodes */
    __forceinline static void swap(Node* a, size_t i, Node* b, size_t j)
    {
      assert(i<4 && j<4);
      std::swap(a->children[i],b->children[j]);
      std::swap(a->lower_x[i],b->lower_x[j]);
      std::swap(a->lower_y[i],b->lower_y[j]);
      std::swap(a->lower_z[i],b->lower_z[j]);
      std::swap(a->upper_x[i],b->upper_x[j]);
      std::swap(a->upper_y[i],b->upper_y[j]);
      std::swap(a->upper_z[i],b->upper_z[j]);
    }

    /*! compacts a node (moves empty children to the end) */
    __forceinline static void compact(Node* a)
    {
      /* find right most filled node */
      ssize_t j=4;
      for (j=j-1; j>=0; j--)
        if (a->child(j) != emptyNode)
          break;

      /* replace empty nodes with filled nodes */
      for (ssize_t i=0; i<j; i++) {
        if (a->child(i) == emptyNode) {
          swap(a,i,a,j);
          for (j=j-1; j>i; j--)
            if (a->child(j) != emptyNode)
              break;
        }
      }
    }
    
  public:

    /*! BVH4 default constructor. */
    BVH4 (const PrimitiveType& primTy, void* geometry = NULL);

    /*! BVH4 destruction */
    ~BVH4 ();

    /*! BVH4 instantiations */
    static Accel* BVH4Triangle1(Scene* scene);
    static Accel* BVH4Triangle4(Scene* scene);
    static Accel* BVH4Triangle8(Scene* scene);
    static Accel* BVH4Triangle1v(Scene* scene);
    static Accel* BVH4Triangle4v(Scene* scene);
    static Accel* BVH4Triangle4i(Scene* scene);
    
    static Accel* BVH4BVH4Triangle1Morton(Scene* scene);
    static Accel* BVH4BVH4Triangle1ObjectSplit(Scene* scene);
    static Accel* BVH4BVH4Triangle4ObjectSplit(Scene* scene);
    static Accel* BVH4BVH4Triangle1vObjectSplit(Scene* scene);
    static Accel* BVH4BVH4Triangle4vObjectSplit(Scene* scene);

    static Accel* BVH4Triangle1SpatialSplit(Scene* scene);
    static Accel* BVH4Triangle4SpatialSplit(Scene* scene);
    static Accel* BVH4Triangle8SpatialSplit(Scene* scene);
    static Accel* BVH4Triangle1ObjectSplit(Scene* scene);
    static Accel* BVH4Triangle4ObjectSplit(Scene* scene);
    static Accel* BVH4Triangle8ObjectSplit(Scene* scene);
    static Accel* BVH4Triangle1vObjectSplit(Scene* scene);
    static Accel* BVH4Triangle4vObjectSplit(Scene* scene);
    static Accel* BVH4Triangle4iObjectSplit(Scene* scene);

    static Accel* BVH4Triangle1ObjectSplit(TriangleMeshScene::TriangleMesh* mesh);
    static Accel* BVH4Triangle4ObjectSplit(TriangleMeshScene::TriangleMesh* mesh);
    static Accel* BVH4Triangle1vObjectSplit(TriangleMeshScene::TriangleMesh* mesh);
    static Accel* BVH4Triangle4vObjectSplit(TriangleMeshScene::TriangleMesh* mesh);
    static Accel* BVH4Triangle4Refit(TriangleMeshScene::TriangleMesh* mesh);

    /*! initializes the acceleration structure */
    void init (size_t numPrimitives = 0);

    /*! Clears the barrier bits of a subtree. */
    void clearBarrier(NodeRef& node);

    Ref<LinearAllocatorPerThread> alloc;

    __forceinline Node* allocNode(size_t thread) {
      Node* node = (Node*) alloc->malloc(thread,sizeof(Node),1 << 7); node->clear(); return node;
    }

    __forceinline char* allocPrimitiveBlocks(size_t thread, size_t num) {
      return (char*) alloc->malloc(thread,num*primTy.bytes,1 << 6);
    }

    /*! Encodes a node */
    __forceinline NodeRef encodeNode(Node* node) { 
      return NodeRef((size_t) node);
    }
    
    /*! Encodes a leaf */
    __forceinline NodeRef encodeLeaf(char* tri, size_t num) {
      assert(!((size_t)tri & align_mask)); 
      return NodeRef((size_t)tri | (1+min(num,(size_t)maxLeafBlocks)));
    }

  public:
    
    /*! calculates the amount of bytes allocated */
    size_t bytesAllocated() 
    {
      if (nodes || primitives)
        return bytesNodes+bytesPrimitives+numVertices*sizeof(Vec3fa);
      else
        return alloc->bytes()+numVertices*sizeof(Vec3fa);
    }

  public:
    const PrimitiveType& primTy;       //!< primitive type stored in the BVH
    void* geometry;                    //!< pointer to additional data for primitive intersector
    NodeRef root;                      //!< Root node
    size_t numPrimitives;
    size_t numVertices;

    /*! data arrays for fast builders */
  public:
    void* nodes;
    size_t bytesNodes;
    void* primitives;
    size_t bytesPrimitives;
    std::vector<BVH4*> objects;
  };

  typedef void (*createTriangleMeshAccelTy)(TriangleMeshScene::TriangleMesh* mesh, BVH4*& accel, Builder*& builder);
  typedef Builder* (*BVH4BuilderTopLevelFunc)(BVH4* accel, Scene* scene, const createTriangleMeshAccelTy createTriangleMeshAccel);

#define DECLARE_TOPLEVEL_BUILDER(symbol)                                         \
  namespace isa   { extern Builder* symbol(BVH4* accel, Scene* scene, const createTriangleMeshAccelTy createTriangleMeshAccel); } \
  namespace sse41 { extern Builder* symbol(BVH4* accel, Scene* scene, const createTriangleMeshAccelTy createTriangleMeshAccel); } \
  namespace avx   { extern Builder* symbol(BVH4* accel, Scene* scene, const createTriangleMeshAccelTy createTriangleMeshAccel); } \
  namespace avx2  { extern Builder* symbol(BVH4* accel, Scene* scene, const createTriangleMeshAccelTy createTriangleMeshAccel); } \
  BVH4BuilderTopLevelFunc symbol;
}

#endif

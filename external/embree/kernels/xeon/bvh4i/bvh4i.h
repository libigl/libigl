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

#include "common/alloc.h"
#include "common/accel.h"
#include "common/scene.h"
#include "geometry/primitive.h"

namespace embree
{
  /*! Multi BVH with 4 children. Each node stores the bounding box of
   * it's 4 children as well as a 4 child indices. */
  class BVH4i : public Bounded
  {
	  ALIGNED_CLASS;
  public:

    /*! forward declaration of node type */
    struct Node;

    /*! branching width of the tree */
    static const size_t N = 4;

    /*! Masks the bits that store the number of items per leaf. */
    static const unsigned offset_mask = 0xFFFFFFFF << 6;
    static const unsigned barrier_mask = 1<<31;
    static const unsigned leaf_mask = 1<<5;  
    static const unsigned items_mask = leaf_mask-1;  
    
    /*! Empty node */
    static const unsigned emptyNode = leaf_mask;

    /*! Invalid node */
    //static const unsigned invalidNode = leaf_mask;
    static const unsigned invalidNode = 0xFFFFFFE0;

    /*! Maximal depth of the BVH. */
    static const size_t maxBuildDepth = 32;
    static const size_t maxBuildDepthLeaf = maxBuildDepth+16;
    static const size_t maxDepth = maxBuildDepthLeaf+maxBuildDepthLeaf+maxBuildDepth; // this depth makes tree rotations of top part of tree safe
    
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
      __forceinline NodeRef (unsigned int id) : id(id) { }

      /*! Cast to unsigned */
      __forceinline operator unsigned int() const { return id; }

      /*! Clears the barrier bit. */
      __forceinline void setBarrier() { id |= barrier_mask; }
      
      /*! Clears the barrier bit. */
      __forceinline void clearBarrier() { id &= ~barrier_mask; }

      /*! Checks if this is an barrier. A barrier tells the top level tree rotation of how deep to enter the tree. */
      __forceinline unsigned int isBarrier() const { return id & barrier_mask; }
     
      /*! checks if this is a leaf */
      __forceinline unsigned int isLeaf() const { return id & leaf_mask; }
      
      /*! checks if this is a node */
      __forceinline unsigned int isNode() const { return (id & leaf_mask) == 0; }
      
      /*! returns node pointer */
      __forceinline       Node* node(      void* base) const { assert(isNode()); return (      Node*)((      char*)base + id); }
      __forceinline const Node* node(const void* base) const { assert(isNode()); return (const Node*)((const char*)base + id); }
      
      /*! returns leaf pointer */
      __forceinline const char* leaf(const void* base, size_t& num) const {
        assert(isLeaf());
        num = id & items_mask;
        return (const char*)base + (id & offset_mask);
      }
      
    private:
      unsigned int id;
    };

    /*! BVH4i Node */
    struct Node
    {
      /*! Clears the node. */
      __forceinline void clear()  {
        lower_x = lower_y = lower_z = pos_inf; 
        upper_x = upper_y = upper_z = neg_inf;
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
        const Vec3fa lower(reduce_min(lower_x),reduce_min(lower_y),reduce_min(lower_z));
        const Vec3fa upper(reduce_max(upper_x),reduce_max(upper_y),reduce_max(upper_z));
        return BBox3f(lower,upper);
      }

      /*! Returns bounds of specified child. */
      __forceinline BBox3f bounds(size_t i) const {
        Vec3fa lower(lower_x[i],lower_y[i],lower_z[i]);
        Vec3fa upper(upper_x[i],upper_y[i],upper_z[i]);
        return BBox3f(lower,upper);
      }

      /*! Returns bounds of all children */
      __forceinline void bounds(BBox<ssef>& bounds0, BBox<ssef>& bounds1, BBox<ssef>& bounds2, BBox<ssef>& bounds3) const {
        transpose(lower_x,lower_y,lower_z,ssef(zero),bounds0.lower,bounds1.lower,bounds2.lower,bounds3.lower);
        transpose(upper_x,upper_y,upper_z,ssef(zero),bounds0.upper,bounds1.upper,bounds2.upper,bounds3.upper);
      }

      /*! Returns reference to specified child */
      __forceinline       NodeRef& child(size_t i)       { return children[i]; }
      __forceinline const NodeRef& child(size_t i) const { return children[i]; }

      /*! Returns number of valid children */
      __forceinline size_t numValidChildren()  {
	size_t valid = 0;
	for (size_t i=0;i<4;i++)
	  if (children[i] != emptyNode)
	    valid++;
	return valid;
      }

    public:
      ssef lower_x;           //!< X dimension of lower bounds of all 4 children.
      ssef upper_x;           //!< X dimension of upper bounds of all 4 children.
      ssef lower_y;           //!< Y dimension of lower bounds of all 4 children.
      ssef upper_y;           //!< Y dimension of upper bounds of all 4 children.
      ssef lower_z;           //!< Z dimension of lower bounds of all 4 children.
      ssef upper_z;           //!< Z dimension of upper bounds of all 4 children.
      NodeRef children[4];    //!< Pointer to the 4 children (can be a node or leaf)
      unsigned int data[4];
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
    
    /*! compacts a node */
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
    BVH4i (const PrimitiveType& primTy, void* geometry = NULL)
      : primTy(primTy), geometry(geometry), maxLeafPrims(maxLeafBlocks*primTy.blockSize), root(emptyNode), qbvh(NULL), accel(NULL)
    {
      alloc_nodes = new LinearAllocatorPerThread;
      alloc_tris  = new LinearAllocatorPerThread;
    }

    /*! BVH4i instantiations */
    static Accel* BVH4iTriangle1(Scene* scene);
    static Accel* BVH4iTriangle4(Scene* scene);
    static Accel* BVH4iTriangle4i(Scene* scene);
    static Accel* BVH4iTriangle1_v1(Scene* scene);
    static Accel* BVH4iTriangle1_v2(Scene* scene);
    static Accel* BVH4iTriangle1_morton(Scene* scene);
    static Accel* BVH4iTriangle1_morton_enhanced(Scene* scene);

    static Accel* BVH4iTriangle8(Scene* scene);

    static Accel* BVH4iTriangle1(TriangleMeshScene::TriangleMesh* mesh);
    static Accel* BVH4iTriangle4(TriangleMeshScene::TriangleMesh* mesh);
    static Accel* BVH4iTriangle1v(TriangleMeshScene::TriangleMesh* mesh);
    static Accel* BVH4iTriangle4v(TriangleMeshScene::TriangleMesh* mesh);
    static Accel* BVH4iTriangle4i(TriangleMeshScene::TriangleMesh* mesh);

    /*! initializes the acceleration structure */
    void init (size_t numNodes = 0, size_t numPrimitives = 0);

    /*! Clears the barrier bits. */
    void clearBarrier(NodeRef& node);

    /*! Calculates the SAH of the BVH */
    float sah ();

  public:

    /*! encodes a node */
    __forceinline NodeRef encodeNode(Node* node) 
    { 
      ssize_t ofs = (size_t)node-(size_t)nodePtr();
      assert(ofs >= 0);
      assert((ofs & ~(size_t)offset_mask) == 0);
      return NodeRef(ofs);
    }
    
    /*! encodes a leaf */
    __forceinline NodeRef encodeLeaf(char* tri, size_t num) 
    {
      ssize_t ofs = (size_t)tri-(size_t)triPtr();
      assert(ofs >= 0);
      assert((ofs & ~(size_t)offset_mask) == 0);
#if defined(_DEBUG)
      if (num > (size_t)maxLeafBlocks) throw std::runtime_error("ERROR: Loosing triangles during build.");
#else
      if (num > (size_t)maxLeafBlocks) std::cout << "WARNING: Loosing triangles during build." << std::endl;
#endif
      return NodeRef(ofs | leaf_mask | min(num,(size_t)maxLeafBlocks));
    }

    /*! Data of the BVH */
  public:
    const size_t maxLeafPrims;          //!< maximal number of triangles per leaf
    NodeRef root;                      //!< Root node (can also be a leaf).

    const PrimitiveType& primTy;   //!< triangle type stored in BVH
    void* geometry;                    //!< pointer to geometry for intersection

 /*! Memory allocation */
  public:
    Ref<LinearAllocatorPerThread> alloc_nodes;
    Ref<LinearAllocatorPerThread> alloc_tris;

    __forceinline Node* allocNode(size_t thread) {
      Node* node = (Node*) alloc_nodes->malloc(thread,sizeof(Node),1 << 7); node->clear(); return node;
    }

    __forceinline char* allocPrimitiveBlocks(size_t thread, size_t num) {
      return (char*) alloc_tris->malloc(thread,num*primTy.bytes,1 << 6);
    }

    __forceinline       void* nodePtr()       { return qbvh; }
    __forceinline const void* nodePtr() const { return qbvh; }

    __forceinline       void* triPtr()       { return accel; }
    __forceinline const void* triPtr() const { return accel; }

    size_t bytes () const {
      return alloc_nodes->bytes() + alloc_tris->bytes();
    }

    // temporaery hack
    void *qbvh;
    void *accel;

  private:
    float sah (NodeRef& node, const BBox3f& bounds);
  };

  __forceinline std::ostream &operator<<(std::ostream &o, const BVH4i::Node &v)
  {
    o << "lower_x " << v.lower_x << std::endl;
    o << "upper_x " << v.upper_x << std::endl;
    
    o << "lower_y " << v.lower_y << std::endl;
    o << "upper_y " << v.upper_y << std::endl;
    
    o << "lower_z " << v.lower_z << std::endl;
    o << "upper_z " << v.upper_z << std::endl;
    
    o << "children " << *(ssei*)v.children << std::endl;
    o << "data     " << *(ssei*)v.data << std::endl;
    
    return o;
  }

  /* ------------------ */
  /* --- Binary BVH --- */
  /* ------------------ */

  extern Vec3fa initQBVHNode[2];

#define BVH_INDEX_SHIFT 6
#define BVH_ITEMS_MASK   (((unsigned int)1 << BVH_INDEX_SHIFT)-1)
#define BVH_LEAF_MASK    ((unsigned int)1 << 31)
#define BVH_OFFSET_MASK  (~(BVH_ITEMS_MASK | BVH_LEAF_MASK))

  template<class T> 
    __forceinline T bvhItemOffset(const T& children) {
    return (children & ~BVH_LEAF_MASK) >> BVH_INDEX_SHIFT;
  }

  template<class T> 
  __forceinline T bvhItems(const T& children) {
    return children & BVH_ITEMS_MASK;
  }
  
  template<class T> 
  __forceinline T bvhChildren(const T& children) {
    return children & BVH_ITEMS_MASK;
  }

  template<class T> 
  __forceinline T bvhChildID(const T& children) {
    return (children & BVH_OFFSET_MASK) >> BVH_INDEX_SHIFT;
  };

  template<class T> 
  __forceinline T bvhLeaf(const T& children) {
    return (children & BVH_LEAF_MASK);
  };

  class __aligned(32) BVHNode : public BBox3f
  {
  public:
    __forceinline unsigned int isLeaf() const {
      return bvhLeaf(lower.a);
    };

    __forceinline int firstChildID() const {
      return bvhChildID(lower.a);
    };
    __forceinline int items() const {
      return bvhItems(lower.a);
    }
    __forceinline unsigned int itemListOfs() const {
      return bvhItemOffset(lower.a);
    }

    __forceinline unsigned int getData() const {
      return upper.a;
    }

    __forceinline void createLeaf(const unsigned int offset,
				  const unsigned int entries,
				  const unsigned int data = 0) 
    {
      assert(entries > 0 && entries <= 4);
      lower.a = (offset << BVH_INDEX_SHIFT) | BVH_LEAF_MASK | entries;
      upper.a = data;
    }

    __forceinline void createNode(const unsigned int index,			  
				  const unsigned short children = 0,
				  const unsigned int items_subtree = 0) {
      assert((index %2) == 0);
      lower.a = (index << BVH_INDEX_SHIFT) | children;
      upper.a = items_subtree;
    }

#ifdef __AVX2__
    __forceinline void operator=(const BVHNode& v) {     
      const avxi b = load8i((int*)&v);
      store8i((int*)this,b);
    };
#endif
  };

  __forceinline std::ostream &operator<<(std::ostream &o, const embree::BVHNode &v)
  {
    if (v.isLeaf())
      {
	o << "LEAF" << " ";
	o << "offset " << v.itemListOfs() << " ";
	o << "items  " << v.items() << " ";
      }
    else
      {
	o << "NODE" << " ";
	o << "firstChildID " << v.firstChildID() << " children " << v.items() << " ";
      }  
    o << "min [" << v.lower <<"] ";
    o << "max [" << v.upper <<"] ";

    return o;
  } 

  /* ---------------- */
  /* --- QUAD BVH --- */
  /* ---------------- */

#define QBVH_INDEX_SHIFT 7
#define QBVH_LEAF_BIT_SHIFT 5
#define QBVH_LEAF_MASK     ((unsigned int)1 << QBVH_LEAF_BIT_SHIFT)
#define QBVH_ITEMS_MASK   (QBVH_LEAF_MASK-1)
#define QBVH_OFFSET_MASK  (~(QBVH_ITEMS_MASK | QBVH_LEAF_MASK))
#define QBVH_TERMINAL_TOKEN QBVH_LEAF_MASK

  class __aligned(64) QBVHNode
  {
  public:
    float min_x[4];
    float max_x[4];
    float min_y[4];
    float max_y[4];
    float min_z[4];
    float max_z[4];
    unsigned int min_d[4];
    unsigned int max_d[4];
  };

  __forceinline unsigned int qbvhItemOffset(const unsigned int children) {
    return children & BVH_OFFSET_MASK; // 6 bits instead of 7
  }

  __forceinline unsigned int qbvhItemOffsetToID(const unsigned int children) {
    return children >> BVH_INDEX_SHIFT; // 6 bits instead of 7
  }

  __forceinline unsigned int qbvhItems(const unsigned int children) {
    return children & QBVH_ITEMS_MASK; // 6 bits instead of 7
  }

  __forceinline unsigned int qbvhChildID(const unsigned int node) {
    return (node & QBVH_OFFSET_MASK) >> QBVH_INDEX_SHIFT;
  };

  __forceinline QBVHNode *qbvhChildPtr(const QBVHNode * __restrict__ const ptr, const unsigned int node) {
    const unsigned int offset = node & QBVH_OFFSET_MASK;
    return (QBVHNode*)((char*)ptr + offset);
  };

  __forceinline QBVHNode *qbvhChildPtrNoMask(const QBVHNode * __restrict__ const ptr, const unsigned int node) {
    return (QBVHNode*)((char*)ptr + (unsigned long)node);
  };

  __forceinline unsigned int qbvhLeaf(const unsigned int node) {
    return (node & QBVH_LEAF_MASK);
  };

  __forceinline unsigned int qbvhLeaf(const unsigned int node, const unsigned int mask) {
    return (node & mask);
  };

  __forceinline unsigned int qbvhChildren(const unsigned int node) {
    return (node & QBVH_ITEMS_MASK);
  };

  template<class T>
    __forceinline T qbvhCreateNode(const T& nodeID, const T& children) {
    return (nodeID << QBVH_INDEX_SHIFT) | children;
  };

}

#endif

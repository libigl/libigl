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

#ifndef __EMBREE_BVH4I_MIC_H__
#define __EMBREE_BVH4I_MIC_H__

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
  public:

    /*! forward declaration of node type */
    struct Node;

    /*! branching width of the tree */
    static const size_t N = 4;

    /*! Masks the bits that store the number of items per leaf. */
    static const unsigned int offset_mask = 0xFFFFFFFF << 6;
    static const unsigned int barrier_mask = 1<<31;
    static const unsigned int leaf_mask = 1<<5;  
    static const unsigned int items_mask = leaf_mask-1;  
    
    /*! Empty node */
    static const unsigned int emptyNode = leaf_mask;

    /*! Invalid node */
    //static const unsigned invalidNode = leaf_mask;
    static const unsigned int invalidNode = 0xFFFFFFE0;

    /*! Maximal depth of the BVH. */
    static const size_t maxBuildDepth = 26;
    static const size_t maxBuildDepthLeaf = maxBuildDepth+6;
    static const size_t maxDepth = maxBuildDepth + maxBuildDepthLeaf;
    
    /*! Cost of one traversal step. */
    static const int travCost = 1;

    static const size_t hybridSIMDUtilSwitchThreshold = 8;

    /*! References a Node or list of Triangles */
    struct NodeRef
    {
      /*! Default constructor */
      __forceinline NodeRef () {}

      /*! Construction from integer */
      __forceinline NodeRef (unsigned id) : _id(id) { }

      /*! Cast to unsigned */
      __forceinline operator unsigned int() const { return _id; }
     
      /*! checks if this is a leaf */
      __forceinline unsigned int isLeaf() const { return _id & leaf_mask; }

      /*! checks if this is a leaf */
      __forceinline unsigned int isLeaf(const unsigned int mask) const { return _id & mask; }
      
      /*! checks if this is a node */
      __forceinline unsigned isNode() const { return (_id & leaf_mask) == 0; }
      
      /*! returns node pointer */
      __forceinline       Node* node(      void* base) const { return (      Node*)((      char*)base + _id); }
      __forceinline const Node* node(const void* base) const { return (const Node*)((const char*)base + _id); }
      
      /*! returns leaf pointer */
      __forceinline const char* leaf(const void* base, unsigned int& num) const {
        assert(isLeaf());
        num = _id & items_mask;
        return (const char*)base + (_id & offset_mask);
      }

      /*! returns leaf pointer */
      __forceinline const char* leaf(const void* base) const {
        assert(isLeaf());
        return (const char*)base + (_id & offset_mask);
      }

      __forceinline unsigned int offset() const {
        return _id & offset_mask;
      }

      __forceinline unsigned int items() const {
        return _id & items_mask;
      }
      
      __forceinline unsigned int &id() { return _id; }
    private:
      unsigned int _id;
    };

    /*! BVH4i Node */

    struct Node
    {
    public:
      struct NodeStruct {
        float x,y,z;           // x,y, and z coordinates of bounds
        NodeRef child;         // encodes 1 is-leaf bit, 25 offset bits, and 6 num-items bits
      } lower[4], upper[4];    // lower and upper bounds of all 4 children

      /*! Returns bounds of specified child. */
      __forceinline BBox3f bounds(size_t i) const {
        Vec3fa l = *(Vec3fa*)&lower[i];
        Vec3fa u = *(Vec3fa*)&upper[i];
        return BBox3f(l,u);
      }

      __forceinline void setInvalid(size_t i)
      {
	lower[i].x = pos_inf;
	lower[i].y = pos_inf;
	lower[i].z = pos_inf;
	lower[i].child = NodeRef(leaf_mask);

	upper[i].x = neg_inf;
	upper[i].y = neg_inf;
	upper[i].z = neg_inf;
	upper[i].child = NodeRef(0);
      }
      /*! Returns reference to specified child */
      __forceinline       NodeRef& child(size_t i)       { return lower[i].child; }
      __forceinline const NodeRef& child(size_t i) const { return lower[i].child; }


      __forceinline std::ostream& operator<<(std::ostream &o)
      {
	for (size_t i=0;i<4;i++)
	  {
	    o << "lower: [" << lower[i].x << "," << lower[i].y << "," << lower[i].z << "] ";
	    o << "upper: [" << upper[i].x << "," << upper[i].y << "," << upper[i].z << "] ";
	    o << std::endl;
	  }
	return o;
      }

    };



  public:

    /*! BVH4 default constructor. */
    BVH4i (const PrimitiveType& primTy, void* geometry = NULL)
      : primTy(primTy), 
      geometry(geometry), 
      root(emptyNode), 
      qbvh(NULL), 
      accel(NULL),
      size_node(0),
      size_accel(0)
    {
    }

    ~BVH4i();

    /*! BVH4i instantiations */
    static Accel* BVH4iTriangle1ObjectSplitBinnedSAH(Scene* scene);
    static Accel* BVH4iTriangle1ObjectSplitMorton(Scene* scene);
    static Accel* BVH4iTriangle1ObjectSplitEnhancedMorton(Scene* scene);
    static Accel* BVH4iVirtualIntersectors(Scene* scene);
    static Accel* BVH4iTriangle1PreSplitsBinnedSAH(Scene* scene);
    static Accel* BVH4iVirtualGeometryBinnedSAH(Scene* scene);

    /*! Calculates the SAH of the BVH */
    float sah ();

    /*! Data of the BVH */
  public:
    //const size_t maxLeafPrims;          //!< maximal number of triangles per leaf
    NodeRef root;                      //!< Root node (can also be a leaf).

    const PrimitiveType& primTy;   //!< triangle type stored in BVH
    void* geometry;                    //!< pointer to geometry for intersection

 /*! Memory allocation */
  public:

    __forceinline       void* nodePtr()       { return qbvh; }
    __forceinline const void* nodePtr() const { return qbvh; }

    __forceinline       void* triPtr()       { return accel; }
    __forceinline const void* triPtr() const { return accel; }


    size_t bytes () const {
      return size_node+size_accel;
    }

    size_t size_node;
    size_t size_accel;

    Node *qbvh;
    void *accel;


    struct Helper { float x,y,z; int a; }; 

    static Helper initQBVHNode[4];

  private:
    float sah (NodeRef& node, BBox3f bounds);
  };


  __forceinline std::ostream &operator<<(std::ostream &o, const BVH4i::Node &v)
    {
      o << std::endl;
      o << "lower: ";
      for (size_t i=0;i<4;i++) o << "[" << v.lower[i].x << "," << v.lower[i].y << "," << v.lower[i].z << "," << v.lower[i].child <<"] ";
      o << std::endl;
      o << "upper: ";
      for (size_t i=0;i<4;i++) o << "[" << v.upper[i].x << "," << v.upper[i].y << "," << v.upper[i].z << "," << v.upper[i].child <<"] ";
      o << std::endl;
      return o;
    } 


  /* ------------------ */
  /* --- Binary BVH --- */
  /* ------------------ */

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

    
    __forceinline void operator=(const BVHNode& v) {     
      const mic_f v_lower = broadcast4to16f((float*)&v.lower);
      const mic_f v_upper = broadcast4to16f((float*)&v.upper);
      store4f((float*)&lower,v_lower);
      store4f((float*)&upper,v_upper);
    };
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

  typedef BVH4i::Node QBVHNode;

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


  __forceinline mic_f initTriangle1(const mic_f &v0,
				    const mic_f &v1,
				    const mic_f &v2,
				    const mic_i &geomID,
				    const mic_i &primID,
				    const mic_i &mask)
  {
    const mic_f e1 = v0 - v1;
    const mic_f e2 = v2 - v0;	     
    const mic_f normal = lcross_xyz(e1,e2);
    const mic_f _v0 = select(0x8888,cast((__m512i)primID),v0);
    const mic_f _v1 = select(0x8888,cast((__m512i)geomID),v1);
    const mic_f _v2 = select(0x8888,cast((__m512i)mask),v2);
    const mic_f _v3 = select(0x8888,mic_f::zero(),normal);
    const mic_f final = lane_shuffle_gather<0>(_v0,_v1,_v2,_v3);
    return final;
  }

  __forceinline void convertToBVH4Layout(BVHNode *__restrict__ const bptr)
  {
    const mic_i box01 = load16i((int*)(bptr + 0));
    const mic_i box23 = load16i((int*)(bptr + 2));

    const mic_i box_min01 = permute<2,0,2,0>(box01);
    const mic_i box_max01 = permute<3,1,3,1>(box01);

    const mic_i box_min23 = permute<2,0,2,0>(box23);
    const mic_i box_max23 = permute<3,1,3,1>(box23);
    const mic_i box_min0123 = select(0x00ff,box_min01,box_min23);
    const mic_i box_max0123 = select(0x00ff,box_max01,box_max23);

    const mic_m min_d_mask = bvhLeaf(box_min0123) != mic_i::zero();
    const mic_i childID    = bvhChildID(box_min0123)>>2;
    const mic_i min_d_node = qbvhCreateNode(childID,mic_i::zero());
    const mic_i min_d_leaf = (box_min0123 ^ BVH_LEAF_MASK) | QBVH_LEAF_MASK;
    const mic_i min_d      = select(min_d_mask,min_d_leaf,min_d_node);
    const mic_i bvh4_min   = select(0x7777,box_min0123,min_d);
    const mic_i bvh4_max   = box_max0123;
    store16i_nt((int*)(bptr + 0),bvh4_min);
    store16i_nt((int*)(bptr + 2),bvh4_max);
  }

}

#endif

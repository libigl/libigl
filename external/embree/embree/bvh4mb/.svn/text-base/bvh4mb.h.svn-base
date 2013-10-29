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

#ifndef __EMBREE_BVH4MB_H__
#define __EMBREE_BVH4MB_H__

#include "../common/accel.h"
#include "../common/alloc.h"
#include "../geometry/triangle.h"

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
      if (num > maxLeafBlocks) std::cout << "WARNING: Loosing triangles during build." << std::endl;
#endif
      return (BaseNode*)((size_t)tri | (1+min(num,size_t(maxLeafBlocks))));
    }
  };

  /*! Multi BVH with 4 children for motion blur. */
  class BVH4MB : public Accel
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

    /*! A base node is either a node or a list of triangles. */
    typedef BaseNode<Node,alignment> Base;

    /*! Maximal depth of the BVH. */
    static const size_t maxDepth = 32;                     
    
    /*! Maximal number of triangle blocks in a leaf. */
    static const size_t maxLeafBlocks = Base::maxLeafBlocks;    

    /*! Cost of one traversal step. */
    static const int travCost = 1;      
   
    /*! BVH4MB Node */
    struct Node
    {
      /*! Clears the node. */
      __forceinline void clear()  {
        lower_x = lower_y = lower_z = ssef(pos_inf);
        upper_x = upper_y = upper_z = ssef(neg_inf);
        lower_dx = lower_dy = lower_dz = ssef(nan); // initialize with NAN and update during refit
        upper_dx = upper_dy = upper_dz = ssef(nan);
        child[0] = child[1] = child[2] = child[3] = (Base*)Base::empty;
      }

      /*! Sets bounding box and ID of child. */
      __forceinline void set(size_t i, const BBox3f& bounds, Base* childID) {
        lower_x[i] = bounds.lower.x; lower_y[i] = bounds.lower.y; lower_z[i] = bounds.lower.z;
        upper_x[i] = bounds.upper.x; upper_y[i] = bounds.upper.y; upper_z[i] = bounds.upper.z;
        child[i] = childID;
      }

      /*! Sets bounding box and ID of child. */
      __forceinline void set(size_t i, const BBox3f& bounds0, const BBox3f& bounds1) 
      {
        lower_x[i] = bounds0.lower.x; lower_y[i] = bounds0.lower.y; lower_z[i] = bounds0.lower.z;
        upper_x[i] = bounds0.upper.x; upper_y[i] = bounds0.upper.y; upper_z[i] = bounds0.upper.z;

        /*! for empty bounds we have to avoid inf-inf=nan */
        if (unlikely(isEmpty(bounds0))) { 
          lower_dx[i] = lower_dy[i] = lower_dz[i] = zero;
          upper_dx[i] = upper_dy[i] = upper_dz[i] = zero;
        } 
        /*! standard case */
        else {
          const Vector3f dlower = bounds1.lower-bounds0.lower;
          const Vector3f dupper = bounds1.upper-bounds0.upper;
          lower_dx[i] = dlower.x; lower_dy[i] = dlower.y; lower_dz[i] = dlower.z;
          upper_dx[i] = dupper.x; upper_dy[i] = dupper.y; upper_dz[i] = dupper.z;
        }
      }

      /*! tests if the node has valid bounds */
      __forceinline bool hasBounds() const {
		  return lower_dx.i[0] != cast_f2i(float(nan));
      }

      /*! Return bounding box for time 0 */
      __forceinline BBox3f bounds0(size_t i) const {
        return BBox3f(Vector3f(lower_x[i],lower_y[i],lower_z[i]),
                      Vector3f(upper_x[i],upper_y[i],upper_z[i]));
      }

      /*! Return bounding box for time 1 */
      __forceinline BBox3f bounds1(size_t i) const {
        return BBox3f(Vector3f(lower_x[i]+lower_dx[i],lower_y[i]+lower_dy[i],lower_z[i]+lower_dz[i]),
                      Vector3f(upper_x[i]+upper_dx[i],upper_y[i]+upper_dy[i],upper_z[i]+upper_dz[i]));
      }

      /*! Returns bounds of node. */
      __forceinline BBox3f bounds() const {
        return BBox3f(Vector3f(reduce_min(min(lower_x,lower_x+lower_dx)),
                            reduce_min(min(lower_y,lower_y+lower_dy)),
                            reduce_min(min(lower_z,lower_z+lower_dz))),
                      Vector3f(reduce_max(max(upper_x,upper_x+upper_dx)),
                            reduce_max(max(upper_y,upper_y+upper_dy)),
                            reduce_max(max(upper_z,upper_z+upper_dz))));
      }

    public:
      ssef_m lower_x, lower_dx;        //!< X dimension of lower bounds of all 4 children.
      ssef_m upper_x, upper_dx;        //!< X dimension of upper bounds of all 4 children.
      ssef_m lower_y, lower_dy;        //!< Y dimension of lower bounds of all 4 children.
      ssef_m upper_y, upper_dy;        //!< Y dimension of upper bounds of all 4 children.
      ssef_m lower_z, lower_dz;        //!< Z dimension of lower bounds of all 4 children.
      ssef_m upper_z, upper_dz;        //!< Z dimension of upper bounds of all 4 children.
      Base* child[4];                //!< Pointer to the 4 children (can be a node or leaf)
    };

    /*! swap the children of two nodes */
    __forceinline static void swap(Node* a, size_t i, Node* b, size_t j)
    {
      std::swap(a->child  [i],b->child  [j]);
      std::swap(a->lower_x[i],b->lower_x[j]);
      std::swap(a->lower_y[i],b->lower_y[j]);
      std::swap(a->lower_z[i],b->lower_z[j]);
      std::swap(a->upper_x[i],b->upper_x[j]);
      std::swap(a->upper_y[i],b->upper_y[j]);
      std::swap(a->upper_z[i],b->upper_z[j]);
      std::swap(a->lower_dx[i],b->lower_dx[j]);
      std::swap(a->lower_dy[i],b->lower_dy[j]);
      std::swap(a->lower_dz[i],b->lower_dz[j]);
      std::swap(a->upper_dx[i],b->upper_dx[j]);
      std::swap(a->upper_dy[i],b->upper_dy[j]);
      std::swap(a->upper_dz[i],b->upper_dz[j]);
    }

  public:

    /*! BVH4MB default constructor. */
    BVH4MB (RTCGeometry* geom, const TriangleType& trity);

    /*! create function for registry */
    static Accel* create (RTCGeometry* geom, const TriangleType& trity) {
      return new BVH4MB(geom,trity);
    }

    /*! clears the acceleration structure */
    void clear ();

    /*! return name of this acceleration structure */
    const std::string name() const { return "bvh4mb."+trity.name; }

    /*! prints statistics */
    void print();

    /*! Rotates tree to improve SAH cost. */
    size_t rotate(Base* node, size_t depth);

    /*! Sort tree to improve shadow ray performance. */
    float sort(RTCGeometry* geom, Base* node, int maxDepth);

    /*! Propagate bounds for time t0 and time t1 up the tree. */
    std::pair<BBox3f,BBox3f> refit(RTCGeometry* geom, Base* node);

    /*! Data of the BVH */
  public:
    const size_t maxLeafTris;          //!< maximal number of triangles per leaf
    AllocatorPerThread alloc;          //!< allocator for nodes and triangles
    Base* root;                        //!< Root node (can also be a leaf).

  private:
    float statistics(Base* node, float area, size_t& depth);
    float bvhSAH;                      //!< SAH cost of the BVH.
    size_t numNodes;                   //!< Number of internal nodes.
    size_t numLeaves;                  //!< Number of leaf nodes.
    size_t numPrimBlocks;              //!< Number of primitive blocks.
    size_t numPrims;                   //!< Number of primitives.
    size_t depth;                      //!< Depth of the tree.
  };
}

#endif

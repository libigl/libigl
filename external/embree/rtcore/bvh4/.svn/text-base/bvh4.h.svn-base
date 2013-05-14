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

#ifndef __EMBREE_BVH4_H__
#define __EMBREE_BVH4_H__

#include "../common/alloc.h"
#include "../common/basenode.h"
#include "../triangle/triangle.h"

namespace embree
{
  /*! Multi BVH with 4 children. Each node stores the bounding box of
   * it's 4 children as well as a 4 child pointers. */
  class BVH4 : public Accel
  {
  public:

    /*! forward declaration of node type */
    struct Node;

    /*! Number of address bits the Node and Triangle are aligned
        to. Maximally 2^alignment-2 many triangle blocks per leaf are
        supported. */
    static const size_t alignment = 4;

    /*! A base node is either a node or a list of triangles. */
    typedef BaseNode<Node,alignment> Base;

    /*! Maximal depth of the BVH. */
    static const size_t maxDepth = 64;
    static const size_t maxLocalDepth = 32;
    
    /*! Maximal number of triangle blocks in a leaf. */
    static const size_t maxLeafBlocks = Base::maxLeafBlocks;    

    /*! Cost of one traversal step. */
    static const int travCost = 1;                      
    
    /*! BVH4 Node */
    struct Node
    {
      /*! Clears the node. */
      __forceinline void clear()  {
        lower_x = lower_y = lower_z = pos_inf;
        upper_x = upper_y = upper_z = neg_inf;
        child[0] = child[1] = child[2] = child[3] = (Base*)Base::empty;
      }

      /*! Sets bounding box and ID of child. */
      __forceinline void set(size_t i, const BBox3f& bounds, Base* childID) {
        lower_x[i] = bounds.lower.x; lower_y[i] = bounds.lower.y; lower_z[i] = bounds.lower.z;
        upper_x[i] = bounds.upper.x; upper_y[i] = bounds.upper.y; upper_z[i] = bounds.upper.z;
        child[i] = childID;
      }

      /*! Sets bounding box and ID of child. */
      __forceinline BBox3f get(size_t i) const {
        return BBox3f(Vec3f(lower_x[i],lower_y[i],lower_z[i]),Vec3f(upper_x[i],upper_y[i],upper_z[i]));
      }

    public:
      ssef lower_x;           //!< X dimension of lower bounds of all 4 children.
      ssef upper_x;           //!< X dimension of upper bounds of all 4 children.
      ssef lower_y;           //!< Y dimension of lower bounds of all 4 children.
      ssef upper_y;           //!< Y dimension of upper bounds of all 4 children.
      ssef lower_z;           //!< Z dimension of lower bounds of all 4 children.
      ssef upper_z;           //!< Z dimension of upper bounds of all 4 children.
      Base* child[4];         //!< Pointer to the 4 children (can be a node or leaf)
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
    }
    
  public:

    /*! BVH4 default constructor. */
    BVH4 (const TriangleType& trity, const std::string& intTy, const Vec3fa* vertices, size_t numVertices, bool freeVertices);

    /*! BVH4 destructor. */
    ~BVH4 ();

    /*! Query interface to the acceleration structure. */
    Ref<RefCount> query(const char* name);

    /*! Print statistics of the BVH. */
    void print(std::ostream& cout);

    /*! Rotates tree to improve SAH cost. */
    size_t rotate(Base* node, size_t depth);

    /*! Sort tree to improve shadow ray performance. */
    float sort(Base* node, int maxDepth);

    /*! Clears the barrier bits. */
    void clearBarrier(Base*& node);

    /*! Data of the BVH */
  public:
    const TriangleType& trity;         //!< triangle type stored in BVH
    const size_t maxLeafTris;          //!< maximal number of triangles per leaf
    AllocatorPerThread alloc;          //!< allocator for nodes and triangles
    Base* root;                        //!< Root node (can also be a leaf).
    const Vec3fa* vertices;            //!< Pointer to vertex array.
    size_t numVertices;                //!< Number of vertices
    bool freeVertices;                 //!< Should we delete the vertex array?

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

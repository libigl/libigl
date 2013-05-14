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

#ifndef __EMBREE_BVH2_H__
#define __EMBREE_BVH2_H__

#include "../common/alloc.h"
#include "../common/basenode.h"
#include "../triangle/triangle.h"

namespace embree
{
  /*! Binary BVH data structure. Each node stores the bounding box of
   * it's 2 children as well as 2 child pointers. */
  class BVH2 : public Accel
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
    static const size_t maxDepth = 128;
    static const size_t maxLocalDepth = 64;
    
    /*! Maximal number of triangle blocks in a leaf. */
    static const size_t maxLeafBlocks = Base::maxLeafBlocks;    

    /*! Cost of one traversal step. */
    static const int travCost = 1;                      
    
    /*! BVH2 Node. The nodes store for each dimension the lower and
     *  upper bounds of the left and right child in one SSE
     *  vector. This allows for a fast swap of the left and right
     *  bounds inside traversal for rays going from "right to left" in
     *  some dimension. */
    struct Node
    {
      /*! Clears the node. */
      __forceinline void clear()  {
        lower_upper_x = ssef(pos_inf,pos_inf,neg_inf,neg_inf);
        child[0] = child[1] = (Base*)Base::empty;
      }

      /*! Sets bounding box and ID of child. */
      __forceinline void set(size_t i, const BBox3f& bounds, Base* childID) {
        lower_upper_x[i+0] = bounds.lower.x;
        lower_upper_y[i+0] = bounds.lower.y;
        lower_upper_z[i+0] = bounds.lower.z;
        lower_upper_x[i+2] = bounds.upper.x;
        lower_upper_y[i+2] = bounds.upper.y;
        lower_upper_z[i+2] = bounds.upper.z;
        child[i] = childID;
      }

      /*! Returns bounds of specified child. */
      __forceinline BBox3f bounds(size_t i) const {
        return BBox3f(Vec3f(lower_upper_x[i+0],lower_upper_y[i+0],lower_upper_z[i+0]),
                      Vec3f(lower_upper_x[i+2],lower_upper_y[i+2],lower_upper_z[i+2]));
      }

    public:
      ssef lower_upper_x;             //!< left_lower_x, right_lower_x, left_upper_x, right_upper_x
      ssef lower_upper_y;             //!< left_lower_y, right_lower_y, left_upper_y, right_upper_y
      ssef lower_upper_z;             //!< left_lower_z, right_lower_z, left_upper_z, right_upper_z
      Base* child[2];                 //!< Pointers to both children.
      char align[16-2*sizeof(void*)]; //!< Padding to one cacheline (64 bytes)
    };

    /*! swap the children of two nodes */
    __forceinline static void swap(Node* a, size_t i, Node* b, size_t j)
    {
      std::swap(a->child        [i  ],b->child        [j  ]);
      std::swap(a->lower_upper_x[i+0],b->lower_upper_x[j+0]);
      std::swap(a->lower_upper_y[i+0],b->lower_upper_y[j+0]);
      std::swap(a->lower_upper_z[i+0],b->lower_upper_z[j+0]);
      std::swap(a->lower_upper_x[i+2],b->lower_upper_x[j+2]);
      std::swap(a->lower_upper_y[i+2],b->lower_upper_y[j+2]);
      std::swap(a->lower_upper_z[i+2],b->lower_upper_z[j+2]);
    }

  public:

    /*! BVH2 default constructor. */
    BVH2 (const TriangleType& trity, const std::string& intTy, const Vec3fa* vertices, size_t numVertices, bool freeVertices);

    /*! BVH2 destructor. */
    ~BVH2 ();

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
    const TriangleType& trity;     //!< triangle type stored in BVH
    const size_t maxLeafTris;      //!< maximal number of triangles per leaf
    AllocatorPerThread alloc;       //!< allocator for nodes and triangles
    Base* root;                    //!< Root node (can also be a leaf).
    const Vec3fa* vertices;        //!< Pointer to vertex array.
    size_t numVertices;            //!< Number of vertices
    bool freeVertices;             //!< Should we delete the vertex array?

  private:
    float statistics(Base* node, float area, size_t& depth);
    float  bvhSAH;                 //!< SAH cost of the BVH.
    size_t numNodes;               //!< Number of internal nodes.
    size_t numLeaves;              //!< Number of leaf nodes.
    size_t numPrimBlocks;          //!< Number of primitive blocks.
    size_t numPrims;               //!< Number of primitives.
    size_t depth;                  //!< Depth of the tree.
  };
}

#endif

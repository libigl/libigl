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

#ifndef __EMBREE_BVH4_STATISTICS_H__
#define __EMBREE_BVH4_STATISTICS_H__

#include "bvh4.h"

namespace embree
{
  class BVH4Statistics 
  {
    typedef BVH4::Node Node;
    typedef BVH4::NodeRef NodeRef;

  public:

    /* Constructor gathers statistics. */
    BVH4Statistics (BVH4* bvh);

    /*! Convert statistics into a string */
    std::string str();

    /*! memory required to store BVH4 */
    size_t bytesUsed();

  private:
    void statistics(NodeRef node, const BBox3f& bounds, size_t& depth);

  private:
    BVH4* bvh;
    float bvhSAH;                      //!< SAH cost of the BVH4.
    float leafSAH;                      //!< SAH cost of the BVH4.
    size_t numNodes;                   //!< Number of internal nodes.
    size_t numLeaves;                  //!< Number of leaf nodes.
    size_t numPrimBlocks;              //!< Number of primitive blocks.
    size_t numPrims;                   //!< Number of primitives.
    size_t depth;                      //!< Depth of the tree.
  };
}

#endif

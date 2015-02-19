// ======================================================================== //
// Copyright 2009-2014 Intel Corporation                                    //
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

#pragma once

#include "bvh8.h"

namespace embree
{
  class BVH8Statistics 
  {
    typedef BVH8::Node Node;
    typedef BVH8::NodeRef NodeRef;

  public:

    /* Constructor gathers statistics. */
    BVH8Statistics (BVH8* bvh);

    /*! Convert statistics into a string */
    std::string str();

    /*! memory required to store BVH8 */
    size_t bytesUsed();

    /*! returns sah cost */
    float sah() const { return bvhSAH; }

  private:
    void statistics(NodeRef node, const BBox3fa& bounds, size_t& depth);

  private:
    BVH8* bvh;
    float bvhSAH;                      //!< SAH cost of the BVH8.
    float leafSAH;                      //!< SAH cost of the BVH8.
    size_t numNodes;                   //!< Number of internal nodes.
    size_t numLeaves;                  //!< Number of leaf nodes.
    size_t numPrimBlocks;              //!< Number of primitive blocks.
    size_t numPrims;                   //!< Number of primitives.
    size_t depth;                      //!< Depth of the tree.
  };
}

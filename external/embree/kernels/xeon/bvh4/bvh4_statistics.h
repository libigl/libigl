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

#include "bvh4.h"

namespace embree
{
  class BVH4Statistics 
  {
    typedef BVH4::Node AlignedNode;
    typedef BVH4::UnalignedNode UnalignedNode;
    typedef BVH4::NodeRef NodeRef;

  public:

    /* Constructor gathers statistics. */
    BVH4Statistics (BVH4* bvh);

    /*! Convert statistics into a string */
    std::string str();

    float sah() const { return bvhSAH; }

    size_t bytesUsed() const;

  private:
    void statistics(NodeRef node, const float A, size_t& depth);

  private:
    BVH4* bvh;
    float bvhSAH;                      //!< SAH cost.
    size_t numAlignedNodes;            //!< Number of aligned internal nodes.
    size_t numUnalignedNodes;          //!< Number of unaligned internal nodes.
    size_t numAlignedNodesMB;            //!< Number of aligned internal nodes.
    size_t numUnalignedNodesMB;          //!< Number of unaligned internal nodes.
    size_t childrenAlignedNodes;       //!< Number of children of aligned nodes
    size_t childrenUnalignedNodes;     //!< Number of children of unaligned internal nodes.
    size_t childrenAlignedNodesMB;       //!< Number of children of aligned nodes
    size_t childrenUnalignedNodesMB;     //!< Number of children of unaligned internal nodes.
    size_t numLeaves;                  //!< Number of leaf nodes.
    size_t numPrims;                   //!< Number of primitives.
    size_t depth;                      //!< Depth of the tree.
    size_t hash;
  };
}

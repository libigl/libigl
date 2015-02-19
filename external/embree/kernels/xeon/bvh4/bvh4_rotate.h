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
  namespace isa 
  {
    /* BVH4 Tree Rotations. */
    class BVH4Rotate
    {
    public:
      typedef BVH4::Node Node;
      typedef BVH4::NodeRef NodeRef;
      
    public:
      static size_t rotate     (BVH4* bvh, NodeRef parentRef, size_t depth = 1);
      
#if defined(__AVX__)
      static void   restructure(NodeRef ref, size_t depth = 1);
#endif
    };
  }
}

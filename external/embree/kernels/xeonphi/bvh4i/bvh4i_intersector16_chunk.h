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

#include "bvh4i.h"
#include "bvh4i_traversal.h"
#include "common/ray16.h" 

namespace embree
{
  namespace isa
  {
    /*! BVH4i Traverser. 16-wide packet traversal implementation for a BVH4i. */
    template<typename LeafIntersector, bool ENABLE_COMPRESSED_BVH4I_NODES>
    class BVH4iIntersector16Chunk
    {
      /* shortcuts for frequently used types */
      typedef typename BVH4i::NodeRef NodeRef;
      typedef typename BVH4i::Node Node;
      
    public:
      static void intersect(mic_i* valid, BVH4i* bvh, Ray16& ray);
      static void occluded (mic_i* valid, BVH4i* bvh, Ray16& ray);
    };

  }
}

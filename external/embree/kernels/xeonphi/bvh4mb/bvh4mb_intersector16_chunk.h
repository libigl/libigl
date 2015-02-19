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

#include "bvh4mb.h"
#include "common/ray16.h" 

namespace embree
{
  namespace isa
  {
    /*! BVH4mb Traverser. Packet traversal implementation for a BVH4mb. */
    template<typename LeafIntersector>
    class BVH4mbIntersector16Chunk
    {
      /* shortcuts for frequently used types */
      typedef typename BVH4i::NodeRef NodeRef;
      typedef typename BVH4i::Node Node;
      
    public:
      static void intersect(mic_i* valid, BVH4mb* bvh, Ray16& ray);
      static void occluded (mic_i* valid, BVH4mb* bvh, Ray16& ray);
    };
  }
}

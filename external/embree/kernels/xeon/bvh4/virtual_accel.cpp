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

#include "virtual_accel.h"

#include "bvh4/bvh4.h"
#include "bvh4i/bvh4i.h"

namespace embree
{
  VirtualAccel::VirtualAccelObjectType virtual_accel_object_type;

  extern Accel::Intersector1 BVH4VirtualIntersector1;
  extern Accel::Intersector4 BVH4VirtualIntersector4Chunk;
  extern Accel::Intersector8 BVH4VirtualIntersector8Chunk;
  extern Accel::Intersector16 BVH4VirtualIntersector16Chunk;

  extern Accel::Intersector1 BVH4iVirtualIntersector1;
  extern Accel::Intersector4 BVH4iVirtualIntersector4Chunk;
  extern Accel::Intersector8 BVH4iVirtualIntersector8Chunk;

  Builder* BVH4BuilderObjectSplit1 (void* bvh, BuildSource* source, void* geometry, const size_t minLeafSize, const size_t maxLeafSize);

  VirtualAccel::VirtualAccel (const std::string& ty, std::vector<AccelSet*>& accels)
    : source(accels)
  {
    if (ty == "bvh4" || ty == "default")
    {
      intersectors.ptr = accel = new BVH4(virtual_accel_object_type);
      intersectors.intersector1 = BVH4VirtualIntersector1;
      intersectors.intersector4 = BVH4VirtualIntersector4Chunk;
      intersectors.intersector8 = BVH4VirtualIntersector8Chunk;
      intersectors.intersector16 = NULL;
      builder = BVH4BuilderObjectSplit1(accel,&source,&accels,1,1);
    }
    else
      throw std::runtime_error("unknown acceleration structure: \"" + ty + "\"");
  }

  VirtualAccel::~VirtualAccel ()
  {
    delete accel;
    delete builder;
  }

  void VirtualAccel::build (size_t threadIndex, size_t threadCount) 
  {
    builder->build(threadIndex,threadCount);
    bounds = accel->bounds;
  }
}

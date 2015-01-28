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

#ifndef __EMBREE_BVHI_BUILDER_MORTON_ENHANCED_MIC_H__
#define __EMBREE_BVHI_BUILDER_MORTON_ENHANCED_MIC_H__

#include "bvh4i_builder_morton.h"

namespace embree
{
  class BVH4iBuilderMortonEnhanced : public BVH4iBuilderMorton
  {
    ALIGNED_CLASS;
  public:

    /*! Constructor. */
    BVH4iBuilderMortonEnhanced (BVH4i* bvh, BuildSource* source, void* geometry);

    /*! creates the builder */
    static Builder* create (void* accel, BuildSource* source, void* geometry) { 
      return new BVH4iBuilderMortonEnhanced((BVH4i*)accel,source,geometry);
    }

    /*! parallel task to iterate over the triangles */
    TASK_RUN_FUNCTION(BVH4iBuilderMortonEnhanced,build_parallel_morton_enhanced);

    /* build function */
    void build(size_t threadIndex, size_t threadCount);

  protected:

    void buildLocalSAHSubTree(BVHNode &parent,
			      BuildRecord &current,
			      PrimRef *__restrict__ local_node,
			      const unsigned int offset_accel);
   

    void extractTopLevelTree(const size_t index,
			     const Vec3fa &root_diag,
			     BVHNode *__restrict__ const local_node,
			     size_t &nodes);

    void buildTopLevelSAHTree(BVHNode &parent,
			      BuildRecord &current,
			      BVHNode *__restrict__ local_node);


  };
}

#endif

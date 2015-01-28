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

#ifndef __EMBREE_BVHI_BUILDER_MORTON_ENHANCED_H__
#define __EMBREE_BVHI_BUILDER_MORTON_ENHANCED_H__

#include "bvh4i_builder_morton.h"

namespace embree
{
  namespace isa
  {
    class BVH4iBuilderMortonEnhanced : public BVH4iBuilderMorton
    {
      ALIGNED_CLASS;
    public:
      
      /*! Constructor. */
      BVH4iBuilderMortonEnhanced (BVH4i* bvh, BuildSource* source, void* geometry, const size_t minLeafSize = 1, const size_t maxLeafSize = inf);
            
      /*! parallel task to iterate over the triangles */
      TASK_RUN_FUNCTION(BVH4iBuilderMortonEnhanced,build_parallel_morton_enhanced);
      
      /* build function */
      void build(size_t threadIndex, size_t threadCount);
      
    protected:
      
      /* static size_t recurseMortonIDWithSAHSubTrees(SmallBuildRecord &current,  */
      /* 						 const size_t mode,  */
      /* 						 const size_t threadID); */
      
      /* static void buildLocalSAHSubTree(SmallBuildRecord &current, */
      /* 				     BVHNode &parent, */
      /* 				     BVHNode *bvh); */
      
      void buildLocalSAHSubTree(BVHNode &parent,
                                BuildRecord &current,
                                PrimRef *__restrict__ local_node,
                                const unsigned int offset_accel);
      
      /* static void thread_recurseSubMortonSAHTrees(void* data, const size_t threadID, const size_t numThreads); */
      
      
      void extractTopLevelTree(const size_t index,
                               const Vec3fa &root_diag,
                               BVHNode *__restrict__ const local_node,
                               size_t &nodes);
      
      void buildTopLevelSAHTree(BVHNode &parent,
                                BuildRecord &current,
                                BVHNode *__restrict__ local_node);
      
      
    };
  }
}

#endif
  

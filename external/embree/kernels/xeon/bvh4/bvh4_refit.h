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

#ifndef __EMBREE_BVH_REFIT_H__
#define __EMBREE_BVH_REFIT_H__

#include "bvh4.h"

namespace embree
{
  namespace isa
  {
    class BVH4Refit : public Builder
    {
      ALIGNED_CLASS;
    public:
      
      /*! Type shortcuts */
      typedef BVH4::Node    Node;
      typedef BVH4::NodeRef NodeRef;
      
    public:
      
      void build(size_t threadIndex, size_t threadCount);
      
      /*! Constructor. */
      BVH4Refit (BVH4* bvh, Builder* builder, TriangleMeshScene::TriangleMesh* mesh);

      ~BVH4Refit();

      TASK_COMPLETE_FUNCTION(BVH4Refit,refit_sequential);
      
      TASK_RUN_FUNCTION(BVH4Refit,task_refit_parallel);
      TASK_COMPLETE_FUNCTION(BVH4Refit,task_refit_complete);
      
    private:
      size_t annotate_tree_sizes(NodeRef& ref);
      void calculate_refit_roots ();
      
      BBox3f leaf_bounds(NodeRef& ref);
      BBox3f node_bounds(NodeRef& ref);
      BBox3f recurse_bottom(NodeRef& ref);
      BBox3f recurse_top(NodeRef& ref);
      
    private:
      //BuildSource* source;           //!< input geometry
      //void* geometry;                //!< input geometry
      TriangleMeshScene::TriangleMesh* mesh;
      
    public:
      const PrimitiveType& primTy;   //!< primitve type stored in BVH
      TaskScheduler::Task task;      //!< parallel refit task
      
    public:
      Builder* builder;
      BVH4* bvh;                      //!< BVH to refit
      std::vector<NodeRef*> roots;    //!< List of equal sized subtrees for bvh refit
    };
  }
}
#endif
  

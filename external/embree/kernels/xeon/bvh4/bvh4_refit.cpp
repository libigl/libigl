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

#include "bvh4_refit.h"
#include "bvh4_statistics.h"
#include "sys/tasklogger.h"

#include <algorithm>

namespace embree
{
  namespace isa
  {
    static const size_t block_size = 1024;
    
    __forceinline bool compare(const BVH4::NodeRef* a, const BVH4::NodeRef* b)
    {
      int sa = *(size_t*)a->node();
      int sb = *(size_t*)b->node();
      return sa < sb;
    }
    
    BVH4Refit::BVH4Refit (BVH4* bvh, Builder* builder, TriangleMesh* mesh, bool listMode)
      : builder(builder), listMode(listMode), mesh(mesh), primTy(bvh->primTy), bvh(bvh), scheduler(&mesh->parent->lockstep_scheduler) 
    {
      needAllThreads = builder->needAllThreads;
    }

    BVH4Refit::~BVH4Refit () {
      delete builder;
    }
    
    void BVH4Refit::build(size_t threadIndex, size_t threadCount) 
    {
      /* build initial BVH */
      if (builder) {
        builder->build(threadIndex,threadCount);
        if (false) { //bvh->numPrimitives > 50000) { // FIXME: reactivate parallel refit
          annotate_tree_sizes(bvh->root);
          calculate_refit_roots();
          needAllThreads = false;
        } else {
          needAllThreads = false;
        }
        delete builder; builder = NULL;
      }
      
      /* refit BVH */
      double t0 = 0.0;
      if (g_verbose >= 2) {
        std::cout << "refitting BVH4 <" << bvh->primTy.name << "> ... " << std::flush;
        t0 = getSeconds();
      }
      
      /* schedule refit tasks */
      size_t numRoots = roots.size();
      if (numRoots <= 1) {
        size_t taskID = TaskLogger::beginTask(threadIndex,"BVH4Refit::sequential",0);
        refit_sequential(threadIndex,threadCount);
        TaskLogger::endTask(threadIndex,taskID);
      }
      else {
        scheduler->dispatchTask(threadIndex,threadCount,_task_refit_parallel,this,numRoots,"BVH4Refit::parallel");
	bvh->bounds = recurse_top(bvh->root);
      }
      
      if (g_verbose >= 2) {
        double t1 = getSeconds();
        std::cout << "[DONE]" << std::endl;
        std::cout << "  dt = " << 1000.0f*(t1-t0) << "ms, perf = " << 1E-6*double(mesh->numTriangles)/(t1-t0) << " Mprim/s" << std::endl;
        std::cout << BVH4Statistics(bvh).str();
      }
    }
    
    size_t BVH4Refit::annotate_tree_sizes(BVH4::NodeRef& ref)
    {
      if (ref.isNode())
      {
        Node* node = ref.node();
        size_t n = 0;
        for (size_t i=0; i<BVH4::N; i++) {
          BVH4::NodeRef& child = node->child(i);
          if (child == BVH4::emptyNode) continue;
          n += annotate_tree_sizes(child); 
        }
        *((size_t*)node) = n;
        return n;
      }
      else
      {
        size_t num; 
        ref.leaf(num);
        return num;
      }
    }
    
    void BVH4Refit::calculate_refit_roots ()
    {
      if (!bvh->root.isNode()) return;
      
      roots.push_back(&bvh->root);
      std::make_heap (roots.begin(), roots.end(), compare);
      
      while (true)
      {
        std::pop_heap(roots.begin(), roots.end(), compare);
        BVH4::NodeRef* node = roots.back();
        roots.pop_back();
        if (*(size_t*)node->node() < block_size) 
          break;
        
        for (size_t i=0; i<BVH4::N; i++) {
          BVH4::NodeRef* child = &node->node()->child(i);
          if (child->isNode()) {
            roots.push_back(child);
            std::push_heap(roots.begin(), roots.end(), compare);
          }
        }
      }
    }
    
    __forceinline BBox3fa BVH4Refit::leaf_bounds(NodeRef& ref)
    {
      size_t num; char* tri = ref.leaf(num);
      if (unlikely(ref == BVH4::emptyNode)) return empty;
      return bvh->primTy.update(tri,listMode ? -1 : num,mesh);
    }
    
    __forceinline BBox3fa BVH4Refit::node_bounds(NodeRef& ref)
    {
      if (ref.isNode())
        return ref.node()->bounds();
      else
        return leaf_bounds(ref);
    }
    
    BBox3fa BVH4Refit::recurse_bottom(NodeRef& ref)
    {
      /* this is a leaf node */
      if (unlikely(ref.isLeaf()))
        return leaf_bounds(ref);
      
      /* recurse if this is an internal node */
      Node* node = ref.node();
      const BBox3fa bounds0 = recurse_bottom(node->child(0));
      const BBox3fa bounds1 = recurse_bottom(node->child(1));
      const BBox3fa bounds2 = recurse_bottom(node->child(2));
      const BBox3fa bounds3 = recurse_bottom(node->child(3));
      
      /* AOS to SOA transform */
      BBox<sse3f> bounds;
      transpose((ssef&)bounds0.lower,(ssef&)bounds1.lower,(ssef&)bounds2.lower,(ssef&)bounds3.lower,bounds.lower.x,bounds.lower.y,bounds.lower.z);
      transpose((ssef&)bounds0.upper,(ssef&)bounds1.upper,(ssef&)bounds2.upper,(ssef&)bounds3.upper,bounds.upper.x,bounds.upper.y,bounds.upper.z);
      
      /* set new bounds */
      node->lower_x = bounds.lower.x;
      node->lower_y = bounds.lower.y;
      node->lower_z = bounds.lower.z;
      node->upper_x = bounds.upper.x;
      node->upper_y = bounds.upper.y;
      node->upper_z = bounds.upper.z;
      
      /* return merged bounds */
      const float lower_x = reduce_min(bounds.lower.x);
      const float lower_y = reduce_min(bounds.lower.y);
      const float lower_z = reduce_min(bounds.lower.z);
      const float upper_x = reduce_max(bounds.upper.x);
      const float upper_y = reduce_max(bounds.upper.y);
      const float upper_z = reduce_max(bounds.upper.z);
      return BBox3fa(Vec3fa(lower_x,lower_y,lower_z),
                    Vec3fa(upper_x,upper_y,upper_z));
    }
    
    BBox3fa BVH4Refit::recurse_top(NodeRef& ref)
    {
      /* stop here if we encounter a barrier */
      if (unlikely(ref.isBarrier())) {
        ref.clearBarrier();
        return node_bounds(ref);
      }
      
      /* this is a leaf node */
      if (unlikely(ref.isLeaf()))
        return leaf_bounds(ref);
      
      /* recurse if this is an internal node */
      Node* node = ref.node();
      const BBox3fa bounds0 = recurse_top(node->child(0));
      const BBox3fa bounds1 = recurse_top(node->child(1));
      const BBox3fa bounds2 = recurse_top(node->child(2));
      const BBox3fa bounds3 = recurse_top(node->child(3));
      
      /* AOS to SOA transform */
      BBox<sse3f> bounds;
      transpose((ssef&)bounds0.lower,(ssef&)bounds1.lower,(ssef&)bounds2.lower,(ssef&)bounds3.lower,bounds.lower.x,bounds.lower.y,bounds.lower.z);
      transpose((ssef&)bounds0.upper,(ssef&)bounds1.upper,(ssef&)bounds2.upper,(ssef&)bounds3.upper,bounds.upper.x,bounds.upper.y,bounds.upper.z);
      
      /* set new bounds */
      node->lower_x = bounds.lower.x;
      node->lower_y = bounds.lower.y;
      node->lower_z = bounds.lower.z;
      node->upper_x = bounds.upper.x;
      node->upper_y = bounds.upper.y;
      node->upper_z = bounds.upper.z;
      
      /* return merged bounds */
      const float lower_x = reduce_min(bounds.lower.x);
      const float lower_y = reduce_min(bounds.lower.y);
      const float lower_z = reduce_min(bounds.lower.z);
      const float upper_x = reduce_max(bounds.upper.x);
      const float upper_y = reduce_max(bounds.upper.y);
      const float upper_z = reduce_max(bounds.upper.z);
      return BBox3fa(Vec3fa(lower_x,lower_y,lower_z),
                    Vec3fa(upper_x,upper_y,upper_z));
    }
    
    void BVH4Refit::task_refit_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount) 
    {
      NodeRef& ref = *roots[taskIndex];
      recurse_bottom(ref);
      ref.setBarrier();
    }
    
    void BVH4Refit::refit_sequential(size_t threadIndex, size_t threadCount) {
      bvh->bounds = recurse_bottom(bvh->root);
    }

    Builder* BVH4Triangle1MeshBuilderFast  (void* bvh, TriangleMesh* mesh, size_t mode);
    Builder* BVH4Triangle4MeshBuilderFast  (void* bvh, TriangleMesh* mesh, size_t mode);
#if defined(__AVX__)
    Builder* BVH4Triangle8MeshBuilderFast  (void* bvh, TriangleMesh* mesh, size_t mode);
#endif
    Builder* BVH4Triangle1vMeshBuilderFast (void* bvh, TriangleMesh* mesh, size_t mode);
    Builder* BVH4Triangle4vMeshBuilderFast (void* bvh, TriangleMesh* mesh, size_t mode);
    Builder* BVH4Triangle4iMeshBuilderFast (void* bvh, TriangleMesh* mesh, size_t mode);

    Builder* BVH4Triangle1MeshRefitFast (void* accel, TriangleMesh* mesh, size_t mode) { return new BVH4Refit((BVH4*)accel,BVH4Triangle1MeshBuilderFast(accel,mesh,mode),mesh,mode); }
    Builder* BVH4Triangle4MeshRefitFast (void* accel, TriangleMesh* mesh, size_t mode) { return new BVH4Refit((BVH4*)accel,BVH4Triangle4MeshBuilderFast(accel,mesh,mode),mesh,mode); }
#if defined(__AVX__)
    Builder* BVH4Triangle8MeshRefitFast (void* accel, TriangleMesh* mesh, size_t mode) { return new BVH4Refit((BVH4*)accel,BVH4Triangle8MeshBuilderFast(accel,mesh,mode),mesh,mode); }
#endif
    Builder* BVH4Triangle1vMeshRefitFast (void* accel, TriangleMesh* mesh, size_t mode) { return new BVH4Refit((BVH4*)accel,BVH4Triangle1vMeshBuilderFast(accel,mesh,mode),mesh,mode); }
    Builder* BVH4Triangle4vMeshRefitFast (void* accel, TriangleMesh* mesh, size_t mode) { return new BVH4Refit((BVH4*)accel,BVH4Triangle4vMeshBuilderFast(accel,mesh,mode),mesh,mode); }
    Builder* BVH4Triangle4iMeshRefitFast (void* accel, TriangleMesh* mesh, size_t mode) { return new BVH4Refit((BVH4*)accel,BVH4Triangle4iMeshBuilderFast(accel,mesh,mode),mesh,mode); }
  }
}


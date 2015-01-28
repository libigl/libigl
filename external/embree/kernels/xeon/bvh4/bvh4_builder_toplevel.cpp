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

#include "bvh4_builder_toplevel.h"
#include "bvh4_statistics.h"

#include "geometry/triangle1.h"
#include "geometry/triangle4.h"
#include "geometry/triangle1v.h"
#include "geometry/triangle4v.h"

#include "sys/tasklogger.h"

namespace embree
{
  namespace isa
  {
#define BUILD_RECORD_SPLIT_THRESHOLD 512
#define THRESHOLD_FOR_SUBTREE_RECURSION 128
#define MIN_OPEN_SIZE 2000

    std::auto_ptr<BVH4BuilderTopLevel::GlobalState> BVH4BuilderTopLevel::g_state(NULL);

    BVH4BuilderTopLevel::BVH4BuilderTopLevel (BVH4* bvh, Scene* scene, const createTriangleMeshAccelTy createTriangleMeshAccel) 
      : bvh(bvh), objects(bvh->objects), scene(scene), createTriangleMeshAccel(createTriangleMeshAccel) {}
    
    BVH4BuilderTopLevel::~BVH4BuilderTopLevel ()
    {
      for (size_t i=0; i<builders.size(); i++) delete builders[i];
      //for (size_t i=0; i<objects.size();  i++) delete objects[i];
    }

    void BVH4BuilderTopLevel::build(size_t threadIndex, size_t threadCount) 
    {
#if 0
      static int frameID = 0;

      if (frameID == 10)
        TaskLogger::start();

      if (frameID == 11) {
        TaskLogger::stop();
        TaskLogger::store("build.fig");
        exit(1);
      }

      frameID++;
#endif

      /* create global state */
      if (!g_state.get()) 
        g_state.reset(new GlobalState(threadCount));

      /* delete some objects */
      size_t N = scene->size();
      for (size_t i=N; i<objects.size(); i++) {
        delete builders[i]; builders[i] = NULL;
        delete objects[i]; objects[i] = NULL;
      }
      
      /* resize object array if scene got larger */
      if (objects.size() < N) {
        objects.resize(N);
        builders.resize(N);
      }
      
      refs.resize(N);
      nextRef = 0;
      
      /* sequential create of acceleration structures */
      for (size_t i=0; i<N; i++) 
        create_object(i);
      
      /* reset bounds of each thread */
      for (size_t i=0; i<threadCount; i++)
        g_state->thread_bounds[i].reset();
      
      /* parallel build of acceleration structures */
      if (N) TaskScheduler::executeTask(threadIndex,threadCount,_task_build_parallel,this,N,"toplevel_build_parallel");
      //for (size_t i=0; i<N; i++) g_state->thread_bounds[threadIndex].extend(build(threadIndex,threadCount,i));
      
      /* perform builds that need all threads */
      for (size_t i=0; i<allThreadBuilds.size(); i++) {
        g_state->thread_bounds[threadIndex].extend(build(threadIndex,threadCount,allThreadBuilds[i]));
      }
      
      allThreadBuilds.clear();
      
      /* build toplevel BVH */
      build_toplevel(threadIndex,threadCount);
    }
    
    void BVH4BuilderTopLevel::build_toplevel(size_t threadIndex, size_t threadCount)
    {
      /* calculate scene bounds */
      Centroid_Scene_AABB bounds; bounds.reset();
      for (size_t i=0; i<threadCount; i++)
        bounds.extend(g_state->thread_bounds[i]);
      
      /* ignore empty scenes */
      //bvh->clear();
      bvh->bounds = bounds.geometry;
      refs.resize(nextRef);
      if (refs.size() == 0) return;
      
      double t0 = 0.0;
      if (g_verbose >= 2) {
        std::cout << "building BVH4<" << bvh->primTy.name << "> with toplevel SAH builder ... " << std::flush;
        t0 = getSeconds();
      }
      
      /* open all large nodes */
#if 0
      open_sequential();
      refs1.resize(refs.size());
#else
      global_dest = refs.size();
      size_t M = max(size_t(2*global_dest),size_t(MIN_OPEN_SIZE));
      refs .resize(M);
      refs1.resize(M);
      barrier.init(threadCount);
      TaskScheduler::executeTask(threadIndex,threadCount,_task_open_parallel,this,threadCount,"toplevel_open_parallel");
      refs.resize(global_dest);
#endif
      bvh->init(refs.size());

      /* start toplevel build */
      BuildRecord task; 
      task.init(bounds,0,refs.size());
      task.parentNode = (size_t)&bvh->root;
      task.depth = 1;
      
      /* initialize thread-local work stacks */
      for (size_t i=0; i<threadCount; i++)
        g_state->thread_workStack[i].reset();
      
      /* push initial build record to global work stack */
      g_state->global_workStack.reset();
      g_state->global_workStack.push_nolock(task);    
      
      /* work in multithreaded toplevel mode until sufficient subtasks got generated */
      while (g_state->global_workStack.size() < 4*threadCount && g_state->global_workStack.size()+BVH4::N <= SIZE_WORK_STACK) 
      {
        BuildRecord br;
        if (!g_state->global_workStack.pop_nolock_largest(br)) break;
        recurseSAH(0,br,BUILD_TOP_LEVEL,threadIndex,threadCount);
      }
      
      /* now process all created subtasks on multiple threads */
      TaskScheduler::executeTask(threadIndex,threadCount,_task_build_subtrees,this,threadCount,"toplevel_build_subtrees");
      
      if (g_verbose >= 2) {
        double t1 = getSeconds();
        std::cout << "[DONE]" << std::endl;
        std::cout << "  dt = " << 1000.0f*(t1-t0) << "ms" << std::endl;
        std::cout << BVH4Statistics(bvh).str();
      }
    }
    
    void BVH4BuilderTopLevel::create_object(size_t objectID)
    {
      TriangleMeshScene::TriangleMesh* mesh = scene->getTriangleMeshSafe(objectID);
      
      /* verify meshes got deleted properly */
      if (mesh == NULL || mesh->numTimeSteps != 1) {
        assert(objectID < objects.size() && objects[objectID] == NULL);
        assert(objectID < builders.size() && builders[objectID] == NULL);
        return;
      }
      
      /* delete BVH and builder for meshes that are scheduled for deletion */
      if (mesh->state == Geometry::ERASING) {
        delete builders[objectID]; builders[objectID] = NULL;
        delete objects[objectID]; objects[objectID] = NULL;
        return;
      }
      
      /* create BVH and builder for new meshes */
      if (objects[objectID] == NULL)
        createTriangleMeshAccel((TriangleMeshScene::TriangleMesh*)mesh,objects[objectID],builders[objectID]);
      
      /* remember meshes that need all threads to get built */
      if (builders[objectID]->needAllThreads) 
        allThreadBuilds.push_back(objectID);
    }
    
    BBox3f BVH4BuilderTopLevel::build (size_t threadIndex, size_t threadCount, size_t objectID)
    {
      /* ignore if no triangle mesh or not enabled */
      TriangleMeshScene::TriangleMesh* mesh = scene->getTriangleMeshSafe(objectID);
      if (mesh == NULL || !mesh->isEnabled() || mesh->numTimeSteps != 1) 
        return empty;
      
      BVH4*    object  = objects [objectID];
      Builder* builder = builders[objectID];
      assert(builder);
      assert(object);
      
      /* build object if it got modified */
      if (mesh->isModified()) {
        builder->build(threadIndex,threadCount);
        mesh->state = Geometry::ENABLED;
      }
      
      /* create build primitive */
      const BBox3f bounds = object->bounds;
      refs[nextRef++] = BuildRef(bounds,object->root);
      return bounds;
    }
    
    void BVH4BuilderTopLevel::task_build_parallel(size_t threadIndex, size_t threadCount, 
                                                  size_t taskIndex, size_t taskCount,
                                                  TaskScheduler::Event* event_i)
    {
      /* ignore meshes that need all threads */
      size_t objectID = taskIndex;
      if (builders[objectID] && builders[objectID]->needAllThreads) 
        return;
      
      /* build all other meshes */
      BBox3f bounds = build(threadIndex,threadCount,objectID);
      if (!bounds.empty()) 
        g_state->thread_bounds[threadIndex].extend(bounds);
    }
    
    void BVH4BuilderTopLevel::open_sequential()
    {
      size_t N = max(2*refs.size(),size_t(MIN_OPEN_SIZE));
      refs.reserve(N);
      
      std::make_heap(refs.begin(),refs.end());
      while (refs.size()+3 <= N)
      {
        std::pop_heap (refs.begin(),refs.end()); 
        BVH4::NodeRef ref = refs.back().node;
        if (ref.isLeaf()) break;
        refs.pop_back();    
        
        BVH4::Node* node = ref.node();
        for (size_t i=0; i<4; i++) {
          if (node->child(i) == BVH4::emptyNode) continue;
          refs.push_back(BuildRef(node->bounds(i),node->child(i)));
          std::push_heap (refs.begin(),refs.end()); 
        }
      }
    }
    
    void BVH4BuilderTopLevel::task_open_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event)
    {
      size_t N = global_dest;
      size_t M = refs1.size();
      const size_t start0 = (threadIndex+0)*N/threadCount;
      const size_t end0   = (threadIndex+1)*N/threadCount;
      const size_t start1 = (threadIndex+0)*M/threadCount;
      const size_t end1   = (threadIndex+1)*M/threadCount;
      assert(end1-start1 >= end0-start0);
      BuildRef* prefs1 = &refs1[0];
      
      /* copy from refs buffer to refs1 buffer */
      for (size_t i=start0, j=start1; i<end0; i++, j++) 
        refs1[j] = refs[i];
      
      /* create max heap in our set of items */
      size_t start = start1;
      size_t end   = start1+end0-start0;
      std::make_heap(&prefs1[start],&prefs1[end]);
      float max_volume = 0.0f;
      
      while (true) 
      {
        barrier.wait(threadIndex,threadCount);
        if (threadIndex == 0) global_max_volume = 0.0f;
        barrier.wait(threadIndex,threadCount);
        
        /* parallel calculation of maximal volume */
        max_volume = 0.0f;
        if (end+3 <= end1)
          for (size_t i=start; i<end; i++)
            max_volume = max(max_volume,prefs1[i].lower.w);
        
        atomic_max_f32(&global_max_volume,max_volume);
        
        barrier.wait(threadIndex,threadCount);
        max_volume = global_max_volume;
        barrier.wait(threadIndex,threadCount);
        
        /* if maximal volume is 0, all threads are finished */
        if (max_volume == 0.0f) break;
                
        /* open all nodes that are considered large in this iteration */
        while (end+3 <= end1)
        {
		  if (end-start == 0) break;
          std::pop_heap(&prefs1[start],&prefs1[end]); 
          BVH4::NodeRef ref = prefs1[end-1].node;
          float vol = prefs1[end-1].lower.w;
		  if (ref.isLeaf() || vol < 0.5f*max_volume) {
			std::push_heap(&prefs1[start],&prefs1[end]); 
			break;
		  }
          end--;
          
          BVH4::Node* node = ref.node();
          for (size_t i=0; i<4; i++) {
            if (node->child(i) == BVH4::emptyNode) continue;
            prefs1[end++] = BuildRef(node->bounds(i),node->child(i));
            std::push_heap(&prefs1[start],&prefs1[end]); 
          }
        }
      }
      
      if (threadIndex == 0) global_dest = 0;
      barrier.wait(threadIndex,threadCount);
      
      /* copy again back to refs array */
      size_t dest = atomic_add(&global_dest,end-start);
      for (size_t i=start, j=dest; i<end; i++, j++) 
        refs[j] = refs1[i];
    }
    
    void BVH4BuilderTopLevel::task_build_subtrees(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event)
    {
      while (true) 
      {
        BuildRecord br;
        if (!g_state->global_workStack.pop_largest(br)) // FIXME: might loose threads during build
        {
          /* global work queue empty => try to steal from neighboring queues */	  
          bool success = false;
          for (size_t i=0; i<threadCount; i++)
          {
            if (g_state->thread_workStack[(threadIndex+i)%threadCount].pop_smallest(br)) {
              success = true;
              break;
            }
          }
          /* found nothing to steal ? */
          if (!success) break; 
        }
        
        /* process local work queue */
        g_state->thread_workStack[threadIndex].push(br);
        while (g_state->thread_workStack[threadIndex].pop_largest(br))
          recurseSAH(0,br,RECURSE,threadIndex,threadCount);
      }
    }
    
    // =======================================================================================================
    // =======================================================================================================
    // =======================================================================================================
    
    void BVH4BuilderTopLevel::split_sequential(BuildRecord& current, BuildRecord& leftChild, BuildRecord& rightChild)
    {
      /* calculate binning function */
      Mapping2<16> mapping(current.bounds);
      
      /* binning of centroids */
      Binner2<16> binner;
      binner.bin(&refs[0],current.begin,current.end,mapping);
      
      /* find best split */
      Split2 split; 
      binner.best(split,mapping);
      
      /* if we cannot find a valid split, enforce an arbitrary split */
      if (unlikely(split.pos == -1)) split_fallback2(&refs[0],current,leftChild,rightChild);
      
      /* partitioning of items */
      else binner.partition(&refs[0], current.begin, current.end, split, mapping, leftChild, rightChild);
    }
    
    void BVH4BuilderTopLevel::split_parallel( BuildRecord &current,
                                              BuildRecord &leftChild,
                                              BuildRecord &rightChild,
                                              const size_t threadID,
                                              const size_t numThreads)
    {
      /* parallel binning of centroids */
      g_state->parallelBinner.bin(current,&refs[0],&refs1[0],threadID,numThreads);
      
      /* find best split */
      Split2 split; 
      g_state->parallelBinner.best(split);
      
      /* if we cannot find a valid split, enforce an arbitrary split */
      if (unlikely(split.pos == -1)) split_fallback2(&refs[0],current,leftChild,rightChild);
      
      /* parallel partitioning of items */
      else g_state->parallelBinner.partition(&refs1[0],&refs[0],split,leftChild,rightChild,threadID,numThreads);
    }
    
    void BVH4BuilderTopLevel::createLeaf(BuildRecord& current, size_t threadIndex, size_t threadCount)
    {
#if defined(DEBUG)
      if (current.depth > BVH4::maxBuildDepthLeaf) 
        throw std::runtime_error("ERROR: depth limit reached");
#endif
      
      /* return empty node */
      if (current.end-current.begin == 0) {
        *(NodeRef*)current.parentNode = BVH4::emptyNode;
        return;
      }
      
      /* return leaf node */
      if (current.end-current.begin == 1) {
        *(NodeRef*)current.parentNode = refs[current.begin].node;
        return;
      }
      
      /* first split level */
      BuildRecord record0, record1;
      split_fallback2(&refs[0],current,record0,record1);
      
      /* second split level */
      BuildRecord children[4];
      split_fallback2(&refs[0],record0,children[0],children[1]);
      split_fallback2(&refs[0],record1,children[2],children[3]);
      
      /* allocate next four nodes */
      BVH4::Node* node = bvh->allocNode(threadIndex);
      *(NodeRef*)current.parentNode = bvh->encodeNode(node);
      
      /* recurse into each child */
      for (size_t i=0; i<4; i++) 
      {
        children[i].parentNode = (size_t)&node->child(i);
        children[i].depth = current.depth+1;
        createLeaf(children[i],threadIndex,threadCount);
        node->set(i,children[i].bounds.geometry);
      }
      BVH4::compact(node); // move empty nodes to the end
    }  
    
    __forceinline void BVH4BuilderTopLevel::split(BuildRecord& current, BuildRecord& left, BuildRecord& right, const size_t mode, const size_t threadID, const size_t numThreads)
    {
      if (mode == BUILD_TOP_LEVEL && current.items() >= BUILD_RECORD_SPLIT_THRESHOLD)
        return split_parallel(current,left,right,threadID,numThreads);		  
      else
        return split_sequential(current,left,right);
    }
    
    __forceinline void BVH4BuilderTopLevel::recurse(size_t depth, BuildRecord& current, const size_t mode, const size_t threadID, const size_t numThreads)
    {
      if (mode == BUILD_TOP_LEVEL) {
        g_state->global_workStack.push_nolock(current);
      }
      else if (current.items() > THRESHOLD_FOR_SUBTREE_RECURSION) {
        if (!g_state->thread_workStack[threadID].push(current))
          recurseSAH(depth,current,RECURSE,threadID,numThreads);
      }
      else
        recurseSAH(depth,current,RECURSE,threadID,numThreads);
    }
    
    void BVH4BuilderTopLevel::recurseSAH(size_t depth, BuildRecord& task, const size_t mode, const size_t threadID, const size_t numThreads)
    {
      /* return leaf node */
      assert(task.end-task.begin > 0);
      if (unlikely(task.end-task.begin == 1)) {
        *(NodeRef*)task.parentNode = refs[task.begin].node;
        return;
      }
      
      /* create leaf node */
      if (unlikely(task.depth >= BVH4::maxBuildDepth)) {
        createLeaf(task,threadID,numThreads);
        return;
      }
      
      /*! initialize task list */
      BuildRecord childTasks[4];
      childTasks[0] = task;
      size_t numChildren = 1;
      
      /*! split until node is full */
      do {
        
        /*! find best child to split */
        float bestArea = inf; 
        ssize_t bestChild = -1;
        for (size_t i=0; i<numChildren; i++) 
        {
          float A = childTasks[i].sceneArea();
          size_t items = childTasks[i].items();
          if (items > 1 && A <= bestArea) { 
            bestChild = i; 
            bestArea = A; 
          }
        }
        if (bestChild == -1) break;
        
        /*! split best child into left and right child */
        __aligned(64) BuildRecord left, right;
        split(childTasks[bestChild],left,right,mode,threadID,numThreads);
        
        /* add new children left and right */
        left.depth = right.depth = task.depth+1;
        childTasks[bestChild] = childTasks[numChildren-1];
        childTasks[numChildren-1] = left;
        childTasks[numChildren+0] = right;
        numChildren++;
        
      } while (numChildren < 4);
      
      /* recurse */
      BVH4::Node* node = bvh->allocNode(threadID);
      for (ssize_t i=numChildren-1; i>=0; i--) {
        childTasks[i].parentNode = (size_t)&node->child(i);
        recurse(depth+1,childTasks[i],mode,threadID,numThreads);
        node->set(i,childTasks[i].bounds.geometry);
      }
      
      *(NodeRef*)task.parentNode = bvh->encodeNode(node);
    }

    Builder* BVH4BuilderTopLevelFast (BVH4* bvh, Scene* scene, const createTriangleMeshAccelTy createTriangleMeshAccel) {
      return new BVH4BuilderTopLevel(bvh,scene,createTriangleMeshAccel);
    }
  }
}

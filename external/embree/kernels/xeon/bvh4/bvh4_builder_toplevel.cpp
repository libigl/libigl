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
#define MIN_OPEN_SIZE 2000

    BVH4BuilderTopLevel::BVH4BuilderTopLevel (BVH4* bvh, Scene* scene, const createTriangleMeshAccelTy createTriangleMeshAccel) 
      : objects(bvh->objects), scene(scene), createTriangleMeshAccel(createTriangleMeshAccel), BVH4TopLevelBuilderFastT(&scene->lockstep_scheduler,bvh) {}
    
    BVH4BuilderTopLevel::~BVH4BuilderTopLevel ()
    {
      for (size_t i=0; i<builders.size(); i++) 
	delete builders[i];
    }

    void BVH4BuilderTopLevel::build(size_t threadIndex, size_t threadCount) 
    {
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
      
      /* parallel build of acceleration structures */
      if (N) scheduler->dispatchTask(threadIndex,threadCount,_task_build_parallel,this,N,"toplevel_build_parallel");
      
      /* perform builds that need all threads */
      for (size_t i=0; i<allThreadBuilds.size(); i++) {
        build(threadIndex,threadCount,allThreadBuilds[i]);
      }
      allThreadBuilds.clear();
      
      /* ignore empty scenes */
      refs.resize(nextRef);
      
      double t0 = 0.0;
      if (g_verbose >= 2) {
	std::cout << "building BVH4<" << bvh->primTy.name << "> with " << TOSTRING(isa) << "::TopLevel SAH builder ... " << std::flush;
        t0 = getSeconds();
      }

      /* open all large nodes */
      open_sequential();

      prims.resize(refs.size());
      for (size_t i=0; i<refs.size(); i++) {
	prims[i] = PrimRef(refs[i].bounds(),(size_t)refs[i].node);
      }
      
      BVH4TopLevelBuilderFastT::build(threadIndex,threadCount,prims.begin(),prims.size());

      if (g_verbose >= 2) {
	std::cout << "[DONE] " << std::endl;
      }
    }

    void BVH4BuilderTopLevel::create_object(size_t objectID)
    {
      TriangleMesh* mesh = scene->getTriangleMeshSafe(objectID);
      
      /* verify meshes got deleted properly */
      if (mesh == NULL || mesh->numTimeSteps != 1) {
        assert(objectID < objects.size () && objects[objectID] == NULL);
        assert(objectID < builders.size() && builders[objectID] == NULL);
        return;
      }
      
      /* delete BVH and builder for meshes that are scheduled for deletion */
      if (mesh->state == Geometry::ERASING) {
        delete builders[objectID]; builders[objectID] = NULL;
        delete objects [objectID]; objects [objectID] = NULL;
        return;
      }
      
      /* create BVH and builder for new meshes */
      if (objects[objectID] == NULL)
        createTriangleMeshAccel((TriangleMesh*)mesh,objects[objectID],builders[objectID]);
      
      /* remember meshes that need all threads to get built */
      if (builders[objectID]->needAllThreads) 
	allThreadBuilds.push_back(objectID);
    }
    
    void BVH4BuilderTopLevel::build (size_t threadIndex, size_t threadCount, size_t objectID)
    {
      /* ignore if no triangle mesh or not enabled */
      TriangleMesh* mesh = scene->getTriangleMeshSafe(objectID);
      if (mesh == NULL || !mesh->isEnabled() || mesh->numTimeSteps != 1) 
        return;
      
      BVH4*    object  = objects [objectID]; assert(object);
      Builder* builder = builders[objectID]; assert(builder);
            
      /* build object if it got modified */
      if (mesh->isModified()) {
        builder->build(threadIndex,threadCount);
        mesh->state = Geometry::ENABLED;
      }
      
      /* create build primitive */
      if (!object->bounds.empty())
	refs[nextRef++] = BuildRef(object->bounds,object->root);
    }
    
    void BVH4BuilderTopLevel::task_build_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount)
    {
      /* ignore meshes that need all threads */
      size_t objectID = taskIndex;
      if (builders[objectID] && builders[objectID]->needAllThreads) 
        return;
      
      /* build all other meshes */
      build(threadIndex,threadCount,objectID);
    }
    
    void BVH4BuilderTopLevel::open_sequential()
    {
      if (refs.size() == 0)
	return;

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
    
    Builder* BVH4BuilderTopLevelFast (BVH4* bvh, Scene* scene, const createTriangleMeshAccelTy createTriangleMeshAccel) {
      return new BVH4BuilderTopLevel(bvh,scene,createTriangleMeshAccel);
    }
  }
}

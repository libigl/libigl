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

#include "scene.h"

#if !defined(__MIC__)
#include "bvh4/bvh4.h"
#include "bvh8/bvh8.h"
#else
#include "xeonphi/bvh4i/bvh4i.h"
#include "xeonphi/bvh4mb/bvh4mb.h"
#include "xeonphi/bvh4hair/bvh4hair.h"
#endif

namespace embree
{
  Scene::Scene (RTCSceneFlags sflags, RTCAlgorithmFlags aflags)
    : flags(sflags), aflags(aflags), numMappedBuffers(0), is_build(false), needTriangles(false), needVertices(false),
      numTriangles(0), numTriangles2(0), 
      numBezierCurves(0), numBezierCurves2(0), 
      numSubdivPatches(0), numSubdivPatches2(0), 
      numUserGeometries1(0), 
      numIntersectionFilters4(0), numIntersectionFilters8(0), numIntersectionFilters16(0),
      commitCounter(0)
  {
#if !defined(__MIC__)
    lockstep_scheduler.taskBarrier.init(TaskScheduler::getNumThreads());
#else
    lockstep_scheduler.taskBarrier.init(MAX_MIC_THREADS);
#endif
    if (g_scene_flags != -1)
      flags = (RTCSceneFlags) g_scene_flags;

    geometries.reserve(128);

#if defined(__MIC__)
    accels.add( BVH4mb::BVH4mbTriangle1ObjectSplitBinnedSAH(this) );
    accels.add( BVH4i::BVH4iVirtualGeometryBinnedSAH(this) );
    accels.add( BVH4Hair::BVH4HairBinnedSAH(this) );
    accels.add( BVH4i::BVH4iSubdivMeshBinnedSAH(this) );

    if (g_verbose >= 1)
      {
	std::cout << "scene flags: static " << isStatic() << " compact = " << isCompact() << " high quality = " << isHighQuality() << " robust = " << isRobust() << std::endl;
      }

    if (g_tri_accel == "default" || g_tri_accel == "bvh4i")   
      {
	if (g_tri_builder == "default") 
	  {
	    if (isStatic())
	      {
		if (g_verbose >= 1) std::cout << "STATIC BUILDER MODE" << std::endl;
		if ( isCompact() )
		  accels.add(BVH4i::BVH4iTriangle1MemoryConservativeBinnedSAH(this));		    
		else if ( isHighQuality() )
		  accels.add(BVH4i::BVH4iTriangle1ObjectSplitBinnedSAH(this));
		else
		  accels.add(BVH4i::BVH4iTriangle1ObjectSplitBinnedSAH(this));
	      }
	    else
	      {
		if (g_verbose >= 1) std::cout << "DYNAMIC BUILDER MODE" << std::endl;
		accels.add(BVH4i::BVH4iTriangle1ObjectSplitMorton(this));
	      }
	  }
	else
	  {
	    if (g_tri_builder == "sah" || g_tri_builder == "bvh4i" || g_tri_builder == "bvh4i.sah") {
	      accels.add(BVH4i::BVH4iTriangle1ObjectSplitBinnedSAH(this));
	    }
	    else if (g_tri_builder == "fast" || g_tri_builder == "morton") {
	      accels.add(BVH4i::BVH4iTriangle1ObjectSplitMorton(this));
	    }
	    else if (g_tri_builder == "fast_enhanced" || g_tri_builder == "morton.enhanced") {
	      accels.add(BVH4i::BVH4iTriangle1ObjectSplitEnhancedMorton(this));
	    }
	    else if (g_tri_builder == "high_quality" || g_tri_builder == "presplits") {
	      accels.add(BVH4i::BVH4iTriangle1PreSplitsBinnedSAH(this));
	    }
	    else if (g_tri_builder == "compact" ||
		     g_tri_builder == "memory_conservative") {
	      accels.add(BVH4i::BVH4iTriangle1MemoryConservativeBinnedSAH(this));
	    }
	    else if (g_tri_builder == "morton64") {
	      accels.add(BVH4i::BVH4iTriangle1ObjectSplitMorton64Bit(this));
	    }

	    else THROW_RUNTIME_ERROR("unknown builder "+g_tri_builder+" for BVH4i<Triangle1>");
	  }
      }
    else THROW_RUNTIME_ERROR("unknown accel "+g_tri_accel);


#else
    createTriangleAccel();
    //accels.add(BVH4::BVH4Triangle1vMB(this));
    accels.add(BVH4::BVH4Triangle4vMB(this));
    accels.add(BVH4::BVH4UserGeometry(this));
    createHairAccel();
    accels.add(BVH4::BVH4OBBBezier1iMB(this,false));
    createSubdivAccel();

#endif
  }

#if !defined(__MIC__)

  void Scene::createTriangleAccel()
  {
    if (g_tri_accel == "default") 
    {
      if (isStatic()) {
        int mode =  2*(int)isCompact() + 1*(int)isRobust(); 
        switch (mode) {
        case /*0b00*/ 0: 
#if defined (__TARGET_AVX__)
          if (has_feature(AVX))
	  {
            if (isHighQuality()) accels.add(BVH8::BVH8Triangle4SpatialSplit(this)); 
            else                 accels.add(BVH8::BVH8Triangle4ObjectSplit(this)); 
          }
          else 
#endif
          {
            if (isHighQuality()) accels.add(BVH4::BVH4Triangle4SpatialSplit(this));
            else                 accels.add(BVH4::BVH4Triangle4ObjectSplit(this)); 
          }
          break;

        case /*0b01*/ 1: accels.add(BVH4::BVH4Triangle4vObjectSplit(this)); break;
        case /*0b10*/ 2: accels.add(BVH4::BVH4Triangle4iObjectSplit(this)); break;
        case /*0b11*/ 3: accels.add(BVH4::BVH4Triangle4iObjectSplit(this)); break;
        }
      } 
      else 
      {
        int mode =  2*(int)isCompact() + 1*(int)isRobust();
        switch (mode) {
        case /*0b00*/ 0: accels.add(BVH4::BVH4BVH4Triangle4ObjectSplit(this)); break;
        case /*0b01*/ 1: accels.add(BVH4::BVH4BVH4Triangle4vObjectSplit(this)); break;
        case /*0b10*/ 2: accels.add(BVH4::BVH4BVH4Triangle4iObjectSplit(this)); break;
        case /*0b11*/ 3: accels.add(BVH4::BVH4BVH4Triangle4iObjectSplit(this)); break;
        }
      }
    }
    else if (g_tri_accel == "bvh4.bvh4.triangle1")    accels.add(BVH4::BVH4BVH4Triangle1ObjectSplit(this));
    else if (g_tri_accel == "bvh4.bvh4.triangle4")    accels.add(BVH4::BVH4BVH4Triangle4ObjectSplit(this));
    else if (g_tri_accel == "bvh4.bvh4.triangle1v")   accels.add(BVH4::BVH4BVH4Triangle1vObjectSplit(this));
    else if (g_tri_accel == "bvh4.bvh4.triangle4v")   accels.add(BVH4::BVH4BVH4Triangle4vObjectSplit(this));
    else if (g_tri_accel == "bvh4.triangle1")         accels.add(BVH4::BVH4Triangle1(this));
    else if (g_tri_accel == "bvh4.triangle4")         accels.add(BVH4::BVH4Triangle4(this));
    else if (g_tri_accel == "bvh4.triangle1v")        accels.add(BVH4::BVH4Triangle1v(this));
    else if (g_tri_accel == "bvh4.triangle4v")        accels.add(BVH4::BVH4Triangle4v(this));
    else if (g_tri_accel == "bvh4.triangle4i")        accels.add(BVH4::BVH4Triangle4i(this));
#if defined (__TARGET_AVX__)
    else if (g_tri_accel == "bvh4.triangle8")         accels.add(BVH4::BVH4Triangle8(this));
    else if (g_tri_accel == "bvh8.triangle4")         accels.add(BVH8::BVH8Triangle4(this));
    else if (g_tri_accel == "bvh8.triangle8")         accels.add(BVH8::BVH8Triangle8(this));
#endif
    else THROW_RUNTIME_ERROR("unknown triangle acceleration structure "+g_tri_accel);
  }

  void Scene::createHairAccel()
  {
    if (g_hair_accel == "default") 
    {
      if (isStatic()) {
        int mode =  2*(int)isCompact() + 1*(int)isRobust(); 
        switch (mode) {
        case /*0b00*/ 0: accels.add(BVH4::BVH4OBBBezier1v(this,isHighQuality())); break;
        case /*0b01*/ 1: accels.add(BVH4::BVH4OBBBezier1v(this,isHighQuality())); break;
        case /*0b10*/ 2: accels.add(BVH4::BVH4OBBBezier1i(this,isHighQuality())); break;
        case /*0b11*/ 3: accels.add(BVH4::BVH4OBBBezier1i(this,isHighQuality())); break;
        }
      } 
      else 
      {
        int mode =  2*(int)isCompact() + 1*(int)isRobust();
        switch (mode) {
	case /*0b00*/ 0: accels.add(BVH4::BVH4Bezier1v(this)); break;
        case /*0b01*/ 1: accels.add(BVH4::BVH4Bezier1v(this)); break;
        case /*0b10*/ 2: accels.add(BVH4::BVH4Bezier1i(this)); break;
        case /*0b11*/ 3: accels.add(BVH4::BVH4Bezier1i(this)); break;
        }
      }   
    }
    else if (g_hair_accel == "bvh4.bezier1v"    ) accels.add(BVH4::BVH4Bezier1v(this));
    else if (g_hair_accel == "bvh4.bezier1i"    ) accels.add(BVH4::BVH4Bezier1i(this));
    else if (g_hair_accel == "bvh4obb.bezier1v" ) accels.add(BVH4::BVH4OBBBezier1v(this,false));
    else if (g_hair_accel == "bvh4obb.bezier1i" ) accels.add(BVH4::BVH4OBBBezier1i(this,false));
    else THROW_RUNTIME_ERROR("unknown hair acceleration structure "+g_hair_accel);
  }

  void Scene::createSubdivAccel()
  {
    if (g_subdiv_accel == "default") 
    {
      if (isIncoherent(flags)) {
        if (isCompact()) accels.add(BVH4::BVH4SubdivGridLazy(this));
        else             accels.add(BVH4::BVH4SubdivGridEager(this));
      }
      else {
        accels.add(BVH4::BVH4SubdivPatch1Cached(this));
      }
    }
    else if (g_subdiv_accel == "bvh4.subdivpatch1"      ) accels.add(BVH4::BVH4SubdivPatch1(this));
    else if (g_subdiv_accel == "bvh4.subdivpatch1cached") accels.add(BVH4::BVH4SubdivPatch1Cached(this));
    else if (g_subdiv_accel == "bvh4.grid.adaptive"     ) accels.add(BVH4::BVH4SubdivGrid(this));
    else if (g_subdiv_accel == "bvh4.grid.eager"        ) accels.add(BVH4::BVH4SubdivGridEager(this));
    else if (g_subdiv_accel == "bvh4.grid.lazy"         ) accels.add(BVH4::BVH4SubdivGridLazy(this));
    else THROW_RUNTIME_ERROR("unknown subdiv accel "+g_subdiv_accel);
  }

#endif

  Scene::~Scene () 
  {
    for (size_t i=0; i<geometries.size(); i++)
      delete geometries[i];
  }

  unsigned Scene::newUserGeometry (size_t items) 
  {
    Geometry* geom = new UserGeometry(this,items);
    return geom->id;
  }
  
  unsigned Scene::newInstance (Scene* scene) {
    Geometry* geom = new Instance(this,scene);
    return geom->id;
  }

  unsigned Scene::newTriangleMesh (RTCGeometryFlags gflags, size_t numTriangles, size_t numVertices, size_t numTimeSteps) 
  {
    if (isStatic() && (gflags != RTC_GEOMETRY_STATIC)) {
      process_error(RTC_INVALID_OPERATION,"static scenes can only contain static geometries");
      return -1;
    }

    if (numTimeSteps == 0 || numTimeSteps > 2) {
      process_error(RTC_INVALID_OPERATION,"only 1 or 2 time steps supported");
      return -1;
    }
    
    Geometry* geom = new TriangleMesh(this,gflags,numTriangles,numVertices,numTimeSteps);
    return geom->id;
  }

  unsigned Scene::newSubdivisionMesh (RTCGeometryFlags gflags, size_t numFaces, size_t numEdges, size_t numVertices, size_t numEdgeCreases, size_t numVertexCreases, size_t numHoles, size_t numTimeSteps) 
  {
    if (isStatic() && (gflags != RTC_GEOMETRY_STATIC)) {
      process_error(RTC_INVALID_OPERATION,"static scenes can only contain static geometries");
      return -1;
    }

    if (numTimeSteps == 0 || numTimeSteps > 2) {
      process_error(RTC_INVALID_OPERATION,"only 1 or 2 time steps supported");
      return -1;
    }
    
    Geometry* geom = new SubdivMesh(this,gflags,numFaces,numEdges,numVertices,numEdgeCreases,numVertexCreases,numHoles,numTimeSteps);
    return geom->id;
  }


  unsigned Scene::newBezierCurves (RTCGeometryFlags gflags, size_t numCurves, size_t numVertices, size_t numTimeSteps) 
  {
    if (isStatic() && (gflags != RTC_GEOMETRY_STATIC)) {
      process_error(RTC_INVALID_OPERATION,"static scenes can only contain static geometries");
      return -1;
    }

    if (numTimeSteps == 0 || numTimeSteps > 2) {
      process_error(RTC_INVALID_OPERATION,"only 1 or 2 time steps supported");
      return -1;
    }
    
    Geometry* geom = new BezierCurves(this,gflags,numCurves,numVertices,numTimeSteps);
    return geom->id;
  }

  unsigned Scene::add(Geometry* geometry) 
  {
    Lock<AtomicMutex> lock(geometriesMutex);

    if (usedIDs.size()) {
      int id = usedIDs.back(); 
      usedIDs.pop_back();
      geometries[id] = geometry;
      return id;
    } else {
      geometries.push_back(geometry);
      return geometries.size()-1;
    }
  }
  
  void Scene::remove(Geometry* geometry) 
  {
    Lock<AtomicMutex> lock(geometriesMutex);
    usedIDs.push_back(geometry->id);
    geometries[geometry->id] = NULL;
    delete geometry;
  }


  void Scene::task_build_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event) 
  {
    LockStepTaskScheduler::Init init(threadIndex,threadCount,&lockstep_scheduler);
    if (threadIndex == 0) accels.build(threadIndex,threadCount);
  }

  void Scene::build (size_t threadIndex, size_t threadCount) 
  {
    /* all user worker threads properly enter and leave the tasking system */
    LockStepTaskScheduler::Init init(threadIndex,threadCount,&lockstep_scheduler);
    if (threadIndex != 0) return;

    /* allow only one build at a time */
    Lock<MutexSys> lock(mutex);

    if (isStatic() && isBuild()) {
      process_error(RTC_INVALID_OPERATION,"static geometries cannot get committed twice");
      return;
    }

    if (!ready()) {
      process_error(RTC_INVALID_OPERATION,"not all buffers are unmapped");
      return;
    }

    /* verify geometry in debug mode  */
#if 0 && defined(DEBUG) // FIXME: enable
    for (size_t i=0; i<geometries.size(); i++) {
      if (geometries[i]) {
        if (!geometries[i]->verify()) {
          process_error(RTC_INVALID_OPERATION,"invalid geometry specified");
          return;
        }
      }
    }
#endif

    /* select fast code path if no intersection filter is present */
    accels.select(numIntersectionFilters4,numIntersectionFilters8,numIntersectionFilters16);

    /* if user provided threads use them */
    if (threadCount)
      accels.build(threadIndex,threadCount);

    /* otherwise use our own threads */
    else
    {
      TaskScheduler::EventSync event;
      new (&task) TaskScheduler::Task(&event,_task_build_parallel,this,TaskScheduler::getNumThreads(),NULL,NULL,"scene_build");
      TaskScheduler::addTask(-1,TaskScheduler::GLOBAL_FRONT,&task);
      event.sync();
    }

    /* make static geometry immutable */
    if (isStatic()) 
    {
      accels.immutable();
      for (size_t i=0; i<geometries.size(); i++)
        if (geometries[i]) geometries[i]->immutable();
    }

    /* delete geometry that is scheduled for delete */
    for (size_t i=0; i<geometries.size(); i++) 
    {
      Geometry* geom = geometries[i];
      if (geom == NULL || geom->state != Geometry::ERASING) continue;
      remove(geom);
    }

    /* update bounds */
    bounds = accels.bounds;
    intersectors = accels.intersectors;
    is_build = true;

    /* enable only algorithms choosen by application */
    if ((aflags & RTC_INTERSECT1) == 0) {
      intersectors.intersector1.intersect = NULL;
      intersectors.intersector1.occluded = NULL;
    }
    if ((aflags & RTC_INTERSECT4) == 0) {
      intersectors.intersector4.intersect = NULL;
      intersectors.intersector4.occluded = NULL;
    }
    if ((aflags & RTC_INTERSECT8) == 0) {
      intersectors.intersector8.intersect = NULL;
      intersectors.intersector8.occluded = NULL;
    }
    if ((aflags & RTC_INTERSECT16) == 0) {
      intersectors.intersector16.intersect = NULL;
      intersectors.intersector16.occluded = NULL;
    }

    if (g_verbose >= 2) {
      std::cout << "created scene intersector" << std::endl;
      accels.print(2);
      std::cout << "selected scene intersector" << std::endl;
      intersectors.print(2);
    }
    
    /* update commit counter */
    commitCounter++;
  }

  void Scene::write(std::ofstream& file)
  {
    int magick = 0x35238765LL;
    file.write((char*)&magick,sizeof(magick));
    int numGroups = size();
    file.write((char*)&numGroups,sizeof(numGroups));
    for (size_t i=0; i<numGroups; i++) {
      if (geometries[i]) geometries[i]->write(file);
      else { int type = -1; file.write((char*)&type,sizeof(type)); }
    }
  }
}

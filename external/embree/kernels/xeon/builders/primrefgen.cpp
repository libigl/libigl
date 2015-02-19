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

#include "primrefgen.h"

namespace embree
{
  namespace isa
  {
    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

    const size_t single_threaded_primrefgen_threshold = 10000;

    void PrimRefListGen::generate(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimRefBlockAlloc<PrimRef>* alloc, const Scene* scene, GeometryTy ty, size_t numTimeSteps, PrimRefList& prims, PrimInfo& pinfo) {
      PrimRefListGen gen(threadIndex,threadCount,scheduler,alloc,scene,ty,numTimeSteps,prims,pinfo);
    }

    PrimRefListGen::PrimRefListGen(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimRefBlockAlloc<PrimRef>* alloc, const Scene* scene, GeometryTy ty, size_t numTimeSteps, PrimRefList& prims, PrimInfo& pinfo)
      : scene(scene), ty(ty), numTimeSteps(numTimeSteps), alloc(alloc), numPrimitives(0), prims_o(prims), pinfo_o(pinfo)
    {
      /*! calculate number of primitives */
      if ((ty & TRIANGLE_MESH) && (numTimeSteps & 1)) numPrimitives += scene->numTriangles;
      if ((ty & TRIANGLE_MESH) && (numTimeSteps & 2)) numPrimitives += scene->numTriangles2;
      if ((ty & SUBDIV_MESH  ) && (numTimeSteps & 1)) numPrimitives += scene->numSubdivPatches;
      if ((ty & SUBDIV_MESH  ) && (numTimeSteps & 2)) numPrimitives += scene->numSubdivPatches2;
      if ((ty & BEZIER_CURVES) && (numTimeSteps & 1)) numPrimitives += scene->numBezierCurves;
      if ((ty & BEZIER_CURVES) && (numTimeSteps & 2)) numPrimitives += scene->numBezierCurves2;
      if ((ty & USER_GEOMETRY)                      ) numPrimitives += scene->numUserGeometries1;
      
      pinfo.reset();
      if (numPrimitives <= single_threaded_primrefgen_threshold) 
	task_gen_parallel(threadIndex,threadCount,0,1);
      else
	scheduler->dispatchTask(threadIndex,threadCount,_task_gen_parallel,this,threadCount,"build::trirefgen");

      assert(pinfo_o.size() <= numPrimitives);
    }
    
    void PrimRefListGen::task_gen_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount) 
    {
      ssize_t start = (taskIndex+0)*numPrimitives/taskCount;
      ssize_t end   = (taskIndex+1)*numPrimitives/taskCount;
      ssize_t cur   = 0;
      
      PrimInfo pinfo(empty);
      PrimRefList::item* block = prims_o.insert(alloc->malloc(threadIndex)); 
      for (size_t i=0; i<scene->size(); i++) 
      {
	const Geometry* geom = scene->get(i);
        if (geom == NULL || !geom->isEnabled()) continue;
	if (!(geom->type & ty)) continue;
	
	switch (geom->type) 
	{
	  /* handle triangle mesh */
	case TRIANGLE_MESH: {
	  const TriangleMesh* mesh = (const TriangleMesh*)geom;
	  if (mesh->numTimeSteps & numTimeSteps) {
	    ssize_t s = max(start-cur,ssize_t(0));
	    ssize_t e = min(end  -cur,ssize_t(mesh->numTriangles));
	    for (ssize_t j=s; j<e; j++) {
	      BBox3fa bounds = empty;
	      if (!mesh->valid(j,&bounds)) continue;
	      const PrimRef prim(bounds,i,j);
	      pinfo.add(prim.bounds(),prim.center2());
	      if (likely(block->insert(prim))) continue; 
	      block = prims_o.insert(alloc->malloc(threadIndex));
	      block->insert(prim);
	    }
	    cur += mesh->numTriangles;
	  }
	  break;
	}

	  /* handle subdiv mesh */
	case SUBDIV_MESH: {
	  const SubdivMesh* mesh = (const SubdivMesh*)geom;
	  if (mesh->numTimeSteps & numTimeSteps) {
	    ssize_t s = max(start-cur,ssize_t(0));
	    ssize_t e = min(end  -cur,ssize_t(mesh->size()));
	    for (ssize_t j=s; j<e; j++) {
	      BBox3fa bounds = empty;
	      if (!mesh->valid(j,&bounds)) continue;
	      const PrimRef prim(bounds,i,j);
	      pinfo.add(prim.bounds(),prim.center2());
	      if (likely(block->insert(prim))) continue; 
	      block = prims_o.insert(alloc->malloc(threadIndex));
	      block->insert(prim);
	    }
	    cur += mesh->size();
	  }
	  break;
	}

	  /* handle bezier curve set */
	case BEZIER_CURVES: {
	  const BezierCurves* set = (const BezierCurves*)geom;
	  if (set->numTimeSteps & numTimeSteps) {
	    ssize_t s = max(start-cur,ssize_t(0));
	    ssize_t e = min(end  -cur,ssize_t(set->numCurves));
	    for (ssize_t j=s; j<e; j++) {
	      BBox3fa bounds = empty;
	      if (!set->valid(j,&bounds)) continue;
	      const PrimRef prim(bounds,i,j);
	      pinfo.add(prim.bounds(),prim.center2());
	      if (likely(block->insert(prim))) continue; 
	      block = prims_o.insert(alloc->malloc(threadIndex));
	      block->insert(prim);
	    }
	    cur += set->numCurves;
	  }
	  break;
	}

	  /* handle user geometry sets */
	case USER_GEOMETRY: {
	  const UserGeometryBase* set = (const UserGeometryBase*)geom;
	  ssize_t s = max(start-cur,ssize_t(0));
	  ssize_t e = min(end  -cur,ssize_t(set->numItems));
	  for (ssize_t j=s; j<e; j++) {
	    BBox3fa bounds = empty;
	    if (!set->valid(j,&bounds)) continue;
	    const PrimRef prim(bounds,i,j);
	    pinfo.add(prim.bounds(),prim.center2());
	    if (likely(block->insert(prim))) continue; 
	    block = prims_o.insert(alloc->malloc(threadIndex));
	    block->insert(prim);
	  }
	  cur += set->numItems;
	  break;
	}
	}
	if (cur >= end) break;  
      }
      pinfo_o.atomic_extend(pinfo);
    }

    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

    template<typename Ty>
    void PrimRefListGenFromGeometry<Ty>::generate(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimRefBlockAlloc<PrimRef>* alloc, const Ty* geom, PrimRefList& prims, PrimInfo& pinfo) {
      PrimRefListGenFromGeometry(threadIndex,threadCount,scheduler,alloc,geom,prims,pinfo);
    }

    template<typename Ty>
    PrimRefListGenFromGeometry<Ty>::PrimRefListGenFromGeometry(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimRefBlockAlloc<PrimRef>* alloc, const Ty* geom, PrimRefList& prims, PrimInfo& pinfo)
      : geom(geom), alloc(alloc), prims_o(prims), pinfo_o(pinfo)
    {
      pinfo_o.reset();
      if (geom->size() <= single_threaded_primrefgen_threshold)
	task_gen_parallel(threadIndex,threadCount,0,1);
      else
	scheduler->dispatchTask(threadIndex,threadCount,_task_gen_parallel,this,threadCount,"build::primrefgen");
      
      assert(pinfo_o.size() <= geom->size());
    }
    
    template<typename Ty>
    void PrimRefListGenFromGeometry<Ty>::task_gen_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount) 
    {
      ssize_t start = (taskIndex+0)*geom->size()/taskCount;
      ssize_t end   = (taskIndex+1)*geom->size()/taskCount;
      ssize_t cur   = 0;
      
      PrimInfo pinfo(empty);
      PrimRefList::item* block = prims_o.insert(alloc->malloc(threadIndex)); 
      
      for (ssize_t j=start; j<end; j++) 
      {
	BBox3fa bounds = empty;
	if (!geom->valid(j,&bounds)) continue;
	const PrimRef prim(bounds,geom->id,j);
	pinfo.add(prim.bounds(),prim.center2());
	if (likely(block->insert(prim))) continue; 
	block = prims_o.insert(alloc->malloc(threadIndex));
	block->insert(prim);
      }
      pinfo_o.atomic_extend(pinfo);
    }
    
    template class PrimRefListGenFromGeometry<TriangleMesh>;
    template class PrimRefListGenFromGeometry<BezierCurves>;
    template class PrimRefListGenFromGeometry<SubdivMesh>;
    template class PrimRefListGenFromGeometry<UserGeometryBase>;

    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

    void PrimRefArrayGen::generate_sequential(size_t threadIndex, size_t threadCount, const Scene* scene, GeometryTy ty, size_t numTimeSteps, PrimRef* prims_o, PrimInfo& pinfo_o) {
      PrimRefArrayGen(threadIndex,threadCount,NULL,scene,ty,numTimeSteps,prims_o,pinfo_o,false);
    }

    void PrimRefArrayGen::generate_parallel(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, const Scene* scene, GeometryTy ty, size_t numTimeSteps, PrimRef* prims_o, PrimInfo& pinfo_o) {
      PrimRefArrayGen(threadIndex,threadCount,scheduler,scene,ty,numTimeSteps,prims_o,pinfo_o,true);
    }

    PrimRefArrayGen::PrimRefArrayGen(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, const Scene* scene, GeometryTy ty, size_t numTimeSteps, PrimRef* prims_o, PrimInfo& pinfo_o, bool parallel)
      : dst(NULL), scene(scene), ty(ty), numTimeSteps(numTimeSteps), numPrimitives(0), prims_o(prims_o), pinfo_o(pinfo_o)
    {
      /*! calculate number of primitives */
      if ((ty & TRIANGLE_MESH) && (numTimeSteps & 1)) numPrimitives += scene->numTriangles;
      if ((ty & TRIANGLE_MESH) && (numTimeSteps & 2)) numPrimitives += scene->numTriangles2;
      if ((ty & SUBDIV_MESH  ) && (numTimeSteps & 1)) numPrimitives += scene->numSubdivPatches;
      if ((ty & SUBDIV_MESH  ) && (numTimeSteps & 2)) numPrimitives += scene->numSubdivPatches2;
      if ((ty & BEZIER_CURVES) && (numTimeSteps & 1)) numPrimitives += scene->numBezierCurves;
      if ((ty & BEZIER_CURVES) && (numTimeSteps & 2)) numPrimitives += scene->numBezierCurves2;
      if ((ty & USER_GEOMETRY)                      ) numPrimitives += scene->numUserGeometries1;

      /*! parallel generation of primref array */
      if (parallel) 
      {
	/* calculate initial destination for each thread */
	dst = new size_t[threadCount];
	for (size_t i=0; i<threadCount; i++)
	  dst[i] = i*numPrimitives/threadCount;

	/* first try to generate primref array */
	pinfo_o.reset();
	scheduler->dispatchTask(task_task_gen_parallel, this, threadIndex, threadCount);
	assert(pinfo_o.size() <= numPrimitives);

	/* calculate new destinations */
	size_t cnt = 0;
	for (size_t i=0; i<threadCount; i++) {
	  size_t n = dst[i]; dst[i] = cnt; cnt += n;
	}

	/* if primitive got filtered out, run again */
	if (cnt < numPrimitives) 
	{
	  pinfo_o.reset();
	  scheduler->dispatchTask(task_task_gen_parallel, this, threadIndex, threadCount);
	  assert(pinfo_o.size() == cnt);
	}

	delete[] dst; dst = NULL;
      }

      /*! sequential generation of primref array */
      else 
      {
	pinfo_o.reset();
	task_gen_parallel(0,1);
	assert(pinfo_o.size() <= numPrimitives);
      }
    }

    void PrimRefArrayGen::task_gen_parallel(size_t taskIndex, size_t taskCount)
    {
      ssize_t start = (taskIndex+0)*numPrimitives/taskCount;
      ssize_t end   = (taskIndex+1)*numPrimitives/taskCount;
      ssize_t dest   = dst ? dst[taskIndex] : start;
      ssize_t cur   = 0;
      
      PrimInfo pinfo(empty);
      for (size_t i=0; i<scene->size(); i++) 
      {
	const Geometry* geom = scene->get(i);
        if (geom == NULL || !geom->isEnabled()) continue;
	if (!(geom->type & ty)) continue;
	
	switch (geom->type) 
	{
	  /* handle triangle mesh */
	case TRIANGLE_MESH: {
	  const TriangleMesh* mesh = (const TriangleMesh*)geom;
	  if (mesh->numTimeSteps & numTimeSteps) {
	    ssize_t s = max(start-cur,ssize_t(0));
	    ssize_t e = min(end  -cur,ssize_t(mesh->numTriangles));
	    for (ssize_t j=s; j<e; j++) {
	      BBox3fa bounds = empty;
	      if (!mesh->valid(j,&bounds)) continue;
	      const PrimRef prim(bounds,i,j);
	      pinfo.add(prim.bounds(),prim.center2());
	      prims_o[dest++] = prim;
	    }
	    cur += mesh->numTriangles;
	  }
	  break;
	}

	  /* handle subdivision meshes */
	case SUBDIV_MESH: {
	  const SubdivMesh* mesh = (const SubdivMesh*)geom;
	  if (mesh->numTimeSteps & numTimeSteps) {
	    ssize_t s = max(start-cur,ssize_t(0));
	    ssize_t e = min(end  -cur,ssize_t(mesh->size()));
	    for (ssize_t j=s; j<e; j++) {
	      BBox3fa bounds = empty;
	      if (!mesh->valid(j,&bounds)) continue;
	      const PrimRef prim(bounds,i,j);
	      pinfo.add(prim.bounds(),prim.center2());
	      prims_o[dest++] = prim;
	    }
	    cur += mesh->size();
	  }
	  break;
	}
	  
	  /* handle bezier curve set */
	case BEZIER_CURVES: {
	  const BezierCurves* set = (const BezierCurves*)geom;
	  if (set->numTimeSteps & numTimeSteps) {
	    ssize_t s = max(start-cur,ssize_t(0));
	    ssize_t e = min(end  -cur,ssize_t(set->numCurves));
	    for (ssize_t j=s; j<e; j++) {
	      BBox3fa bounds = empty;
	      if (!set->valid(j,&bounds)) continue;
	      const PrimRef prim(bounds,i,j);
	      pinfo.add(prim.bounds(),prim.center2());		
	      prims_o[dest++] = prim;
	    }
	    cur += set->numCurves;
	  }
	  break;
	}
	  
	  /* handle user geometry sets */
	case USER_GEOMETRY: {
	  const UserGeometryBase* set = (const UserGeometryBase*)geom;
	  ssize_t s = max(start-cur,ssize_t(0));
	  ssize_t e = min(end  -cur,ssize_t(set->numItems));
	  for (ssize_t j=s; j<e; j++) {
	    BBox3fa bounds = empty;
	    if (!set->valid(j,&bounds)) continue;
	    const PrimRef prim(bounds,i,j);
	    if (!inFloatRange(prim.bounds())) continue;
	    pinfo.add(prim.bounds(),prim.center2());
	    prims_o[dest++] = prim;
	  }
	  cur += set->numItems;
	  break;
	}
	}
	if (cur >= end) break;  
      }
      pinfo_o.atomic_extend(pinfo);
      if (dst) dst[taskIndex] = dest - dst[taskIndex];
    }

    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

    template<typename Ty>
      void PrimRefArrayGenFromGeometry<Ty>::generate_sequential(size_t threadIndex, size_t threadCount, const Ty* geom, PrimRef* prims_o, PrimInfo& pinfo_o) 
    {
      pinfo_o.reset();
      PrimRefArrayGenFromGeometry gen(geom,prims_o,pinfo_o);
      gen.task_gen_parallel(0,1);
      assert(pinfo_o.size() <= geom->size());
    }
    
    template<typename Ty>
      void PrimRefArrayGenFromGeometry<Ty>::generate_parallel(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, const Ty* geom, PrimRef* prims_o, PrimInfo& pinfo_o) 
    {
      PrimRefArrayGenFromGeometry gen(geom,prims_o,pinfo_o);

      /* calculate initial destination for each thread */
      size_t numPrimitives = geom->size();
      gen.dst = new size_t[threadCount];
      for (size_t i=0; i<threadCount; i++)
	gen.dst[i] = i*numPrimitives/threadCount;
      
      /* first try to generate primref array */
      pinfo_o.reset();
      scheduler->dispatchTask(task_task_gen_parallel, &gen, threadIndex, threadCount);
      assert(pinfo_o.size() <= numPrimitives);
      
      /* calculate new destinations */
      size_t cnt = 0;
      for (size_t i=0; i<threadCount; i++) {
	size_t n = gen.dst[i]; gen.dst[i] = cnt; cnt += n;
      }
      
      /* if primitive got filtered out, run again */
      if (cnt < numPrimitives) 
      {
	pinfo_o.reset();
	scheduler->dispatchTask(task_task_gen_parallel, &gen, threadIndex, threadCount);
	assert(pinfo_o.size() == cnt);
      }
      
      delete[] gen.dst; gen.dst = NULL;
    }

    template<typename Ty>
      PrimRefArrayGenFromGeometry<Ty>::PrimRefArrayGenFromGeometry(const Ty* geom, PrimRef* prims_o, PrimInfo& pinfo_o)
	: geom(geom), prims_o(prims_o), pinfo_o(pinfo_o), dst(NULL) {}

    template<typename Ty>
      void PrimRefArrayGenFromGeometry<Ty>::task_gen_parallel(size_t taskIndex, size_t taskCount)
    {
      ssize_t start = (taskIndex+0)*geom->size()/taskCount;
      ssize_t end   = (taskIndex+1)*geom->size()/taskCount;
      ssize_t dest  = dst ? dst[taskIndex] : start;
      ssize_t cur   = 0;
      
      PrimInfo pinfo(empty);
      for (ssize_t j=start; j<end; j++) {
	BBox3fa bounds = empty;
	if (!geom->valid(j,&bounds)) continue;
	const PrimRef prim(bounds,geom->id,j);
	pinfo.add(prim.bounds(),prim.center2());
	prims_o[dest++] = prim;
      }
      pinfo_o.atomic_extend(pinfo);
      if (dst) dst[taskIndex] = dest - dst[taskIndex];
    }

    template class PrimRefArrayGenFromGeometry<TriangleMesh>;
    template class PrimRefArrayGenFromGeometry<BezierCurves>;
    template class PrimRefArrayGenFromGeometry<SubdivMesh>;
    template class PrimRefArrayGenFromGeometry<UserGeometryBase>;
  }
}


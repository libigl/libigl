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

#include "bezierrefgen.h"

namespace embree
{
  namespace isa
  {
    BezierRefGen::BezierRefGen(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimRefBlockAlloc<BezierPrim>* alloc, const Scene* scene, const size_t numTimeSteps)
      : scene(scene), numTimeSteps(numTimeSteps), alloc(alloc), pinfo(empty)
    {
      /*! parallel stage */
      size_t numTasks = min(threadCount,maxTasks);
      scheduler->dispatchTask(threadIndex,threadCount,_task_gen_parallel,this,numTasks,"build::bezierrefgen");
      
      /*! reduction stage */
      for (size_t i=0; i<numTasks; i++)
	pinfo.merge(pinfos[i]);
    }
    
    void BezierRefGen::task_gen_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount) 
    {
      ssize_t numBezierCurves = numTimeSteps == 2 ? scene->numBezierCurves2 : scene->numBezierCurves;
      ssize_t start = (taskIndex+0)*numBezierCurves/taskCount;
      ssize_t end   = (taskIndex+1)*numBezierCurves/taskCount;
      ssize_t cur   = 0;
      
      PrimInfo pinfo(empty);
      BezierRefList::item* block = prims.insert(alloc->malloc(threadIndex)); 
      for (size_t i=0; i<scene->size(); i++) 
      {
	BezierCurves* geom = (BezierCurves*) scene->get(i);
        if (geom == NULL) continue;
	if (geom->type != BEZIER_CURVES || !geom->isEnabled() || geom->numTimeSteps != numTimeSteps) continue;
	ssize_t gstart = 0;
	ssize_t gend = geom->numCurves;
	ssize_t s = max(start-cur,gstart);
	ssize_t e = min(end  -cur,gend  );
	for (size_t j=s; j<e; j++) 
	{
	  if (!geom->valid(j)) continue;

	  const int ofs = geom->curve(j);
	  Vec3fa p0 = geom->vertex(ofs+0,0);
	  Vec3fa p1 = geom->vertex(ofs+1,0);
	  Vec3fa p2 = geom->vertex(ofs+2,0);
	  Vec3fa p3 = geom->vertex(ofs+3,0);
	  if (numTimeSteps == 2) {
	    p0 = 0.5f*(p0+geom->vertex(ofs+0,1));
	    p1 = 0.5f*(p1+geom->vertex(ofs+1,1));
	    p2 = 0.5f*(p2+geom->vertex(ofs+2,1));
	    p3 = 0.5f*(p3+geom->vertex(ofs+3,1));
	  }
	  const BezierPrim bezier(p0,p1,p2,p3,0,1,i,j,false);
	  pinfo.add(bezier.bounds(),bezier.center());
	  if (likely(block->insert(bezier))) continue; 
	  block = prims.insert(alloc->malloc(threadIndex));
	  block->insert(bezier);
	}
	cur += geom->numCurves;
	if (cur >= end) break;  
      }
      pinfos[taskIndex] = pinfo;
    }
  }
}

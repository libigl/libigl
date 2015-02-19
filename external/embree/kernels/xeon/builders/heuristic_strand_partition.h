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

#pragma once

#include "geometry/bezier1v.h"
#include "builders/primrefalloc.h"
#include "heuristic_fallback.h"

namespace embree
{
  namespace isa
  {
    /*! Tries to split hair into two differently aligned hair strands */
    struct StrandSplit
    {
      struct Split;
      typedef atomic_set<PrimRefBlockT<BezierPrim> > BezierRefList;
      
    public:
      StrandSplit () {}
      
      /*! finds the two hair strands */
      template<bool Parallel>
	static const Split find(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, BezierRefList& curves);
      
    private:
      
      /*! number of tasks */
      static const size_t maxTasks = 32;
      
    public:
      
      /*! stores all information to perform some split */
      struct Split
      {    
	/*! construct an invalid split by default */
	__forceinline Split()
	  : sah(inf), axis0(zero), axis1(zero) {}
	
	/*! constructs specified split */
	__forceinline Split(const float sah, const Vec3fa& axis0, const Vec3fa& axis1)
	  : sah(sah), axis0(axis0), axis1(axis1) {}
	
	/*! calculates standard surface area heuristic for the split */
	__forceinline float splitSAH() const { return sah; }
	
	/*! splitting into two sets */
	template<bool Parallel>
	  void split(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimRefBlockAlloc<BezierPrim>& alloc, 
		     BezierRefList& prims, 
		     BezierRefList& lprims_o, PrimInfo& linfo_o, 
		     BezierRefList& rprims_o, PrimInfo& rinfo_o) const;
	
      private:
	float sah;             //!< SAH cost of the split
	Vec3fa axis0, axis1;   //!< axis the two strands are aligned into
      };
      
    private:
      
      /*! parallel find of split */
      struct TaskFindParallel
      {
	/*! construction executes the task */
	TaskFindParallel(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, BezierRefList& prims);
	
      private:
	
	/*! parallel task function */
	TASK_SET_FUNCTION(TaskFindParallel,task_find_parallel);
	
	/*! parallel bounding calculations */
	TASK_SET_FUNCTION(TaskFindParallel,task_bound_parallel);
	
	/*! state for find stage */
      private:
	BezierRefList::iterator iter0; //!< iterator for find stage
	Vec3fa axis0;
	Vec3fa axis1;
	float task_cos[maxTasks];
	Vec3fa task_axis1[maxTasks];
	
	/*! state for bounding stage */
      private:
	BezierRefList::iterator iter1; //!< iterator for bounding stage 
	size_t task_lnum[maxTasks];
	size_t task_rnum[maxTasks];
	BBox3fa task_lbounds[maxTasks];
	BBox3fa task_rbounds[maxTasks];
	
      public:
	Split split; //!< best split
      };
      
      /*! task for parallel splitting */
      struct TaskSplitParallel
      {
	/*! construction executes the task */
	TaskSplitParallel(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, const Split* split, PrimRefBlockAlloc<BezierPrim>& alloc, 
			  BezierRefList& prims, 
			  BezierRefList& lprims_o, PrimInfo& linfo_o, 
			  BezierRefList& rprims_o, PrimInfo& rinfo_o);
	
      private:
	
	/*! parallel split task function */
	TASK_SET_FUNCTION(TaskSplitParallel,task_split_parallel);
	
	/*! input data */
      private:
	const Split* split;
	PrimRefBlockAlloc<BezierPrim>& alloc;
	BezierRefList prims;
	PrimInfo linfos[maxTasks];
	PrimInfo rinfos[maxTasks];
	
	/*! output data */
      private:
	BezierRefList& lprims_o; 
	PrimInfo& linfo_o;
	BezierRefList& rprims_o;
	PrimInfo& rinfo_o;
      };
    };
  }
}

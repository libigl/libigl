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
#include "common/scene.h"

namespace embree
{
  namespace isa
  {
    struct SpatialSplit
    {
      struct Split;
      typedef atomic_set<PrimRefBlockT<PrimRef> > TriRefList;    //!< list of triangles
      typedef atomic_set<PrimRefBlockT<BezierPrim> > BezierRefList; //!< list of bezier primitives
      
    public:
      
      /*! finds the best split */
      template<bool Parallel>
      static const Split find(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, Scene* scene, BezierRefList& curves, const PrimInfo& pinfo, const size_t logBlockSize);
      
      /*! finds the best split */
      template<bool Parallel>
      static const Split find(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, Scene* scene, TriRefList& curves, const PrimInfo& pinfo, const size_t logBlockSize);
      
    private:
      
      /*! number of bins */
      static const size_t BINS = 16;
      
      /*! number of tasks */
      static const size_t maxTasks = 32;
            
      /*! mapping into bins */
      struct Mapping
      {
      public:
	__forceinline Mapping() {}
	
	/*! calculates the mapping */
	__forceinline Mapping(const PrimInfo& pinfo);
	
	/*! slower but safe binning */
	__forceinline ssei bin(const Vec3fa& p) const;
	
	/*! calculates left spatial position of bin */
	__forceinline float pos(const int bin, const int dim) const;
	
	/*! returns true if the mapping is invalid in some dimension */
	__forceinline bool invalid(const int dim) const;
	
      private:
	ssef ofs,scale;  //!< linear function that maps to bin ID
      };
      
    public:
      
      /*! stores all information required to perform some split */
      struct Split
      {
	/*! construct an invalid split by default */
	__forceinline Split() 
	  : sah(inf), dim(-1), pos(0) {}
	
	/*! constructs specified split */
	__forceinline Split(float sah, int dim, int pos, const Mapping& mapping)
	  : sah(sah), dim(dim), pos(pos), mapping(mapping) {}

	/*! tests if this split is valid */
	__forceinline bool valid() const { return dim != -1; }

	/*! calculates surface area heuristic for performing the split */
	__forceinline float splitSAH() const { return sah; }
	
	/*! splitting into two sets */
	template<bool Parallel>
	  void split(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, 
		     PrimRefBlockAlloc<BezierPrim>& alloc, 
		     Scene* scene, BezierRefList& curves, 
		     BezierRefList& lprims_o, PrimInfo& linfo_o, 
		     BezierRefList& rprims_o, PrimInfo& rinfo_o) const;

	/*! splitting into two sets */
	template<bool Parallel>
	  void split(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, 
		     PrimRefBlockAlloc<PrimRef>& alloc, 
		     Scene* scene, TriRefList& curves, 
		     TriRefList& lprims_o, PrimInfo& linfo_o, 
		     TriRefList& rprims_o, PrimInfo& rinfo_o) const;

	/*! stream output */
	friend std::ostream& operator<<(std::ostream& cout, const Split& split) {
	  return cout << "Split { sah = " << split.sah << ", dim = " << split.dim << ", pos = " << split.pos << "}";
	}
	
      public:
	float sah;          //!< SAH cost of the split
	int   dim;          //!< split dimension
	int   pos;          //!< split position
	Mapping mapping;    //!< mapping into bins
      };
      
    private:
      
      /*! stores all binning information */
      struct __aligned(64) BinInfo
      {
	BinInfo();
	
	/*! bins an array of bezier primitives */
	void bin (Scene* scene, const BezierPrim* prims, size_t N, const PrimInfo& pinfo, const Mapping& mapping);

	/*! bins an array of triangles */
	void bin(Scene* scene, const PrimRef* prims, size_t N, const PrimInfo& pinfo, const Mapping& mapping);
	
	/*! bins a list of bezier primitives */
	void bin(Scene* scene, BezierRefList& prims, const PrimInfo& pinfo, const Mapping& mapping);
	
	/*! bins a list of triangles */
	void bin(Scene* scene, TriRefList& prims, const PrimInfo& pinfo, const Mapping& mapping);
	
	/*! merges in other binning information */
	void merge (const BinInfo& other);
	
	/*! finds the best split by scanning binning information */
	Split best(const PrimInfo& pinfo, const Mapping& mapping, const size_t logBlockSize);
	
      private:
	BBox3fa bounds[BINS][4];  //!< geometry bounds for each bin in each dimension
	ssei    numBegin[BINS];   //!< number of primitives starting in bin
	ssei    numEnd[BINS];     //!< number of primitives ending in bin
      };
      
      /*! task for parallel binning */
      template<typename List>
      struct TaskBinParallel
      {
	/*! construction executes the task */
	TaskBinParallel(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, Scene* scene, List& prims, const PrimInfo& pinfo, const Mapping& mapping, const size_t logBlockSize);
	
      private:
	
	/*! parallel binning */
	TASK_SET_FUNCTION(TaskBinParallel,task_bin_parallel);
	
      private:
	Scene* scene;
	typename List::iterator iter; 
	PrimInfo pinfo;
	Mapping mapping;
	BinInfo binners[maxTasks];
	
      public:
	Split split; //!< best split
      };
      
      /*! task for parallel splitting */ 
      template<typename Prim>
      struct TaskSplitParallel
      {
	typedef atomic_set<PrimRefBlockT<Prim> > List;

	/*! construction executes the task */
	TaskSplitParallel(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, const Split* split, PrimRefBlockAlloc<Prim>& alloc, 
			  Scene* scene, List& prims, 
			  List& lprims_o, PrimInfo& linfo_o, 
			  List& rprims_o, PrimInfo& rinfo_o);
	
      private:
	
	/*! parallel split task function */
	TASK_SET_FUNCTION(TaskSplitParallel,task_split_parallel);
	
	/*! input data */
      private:
	const Split* split;
	PrimRefBlockAlloc<Prim>& alloc;
	Scene* scene;
	List prims;
	PrimInfo linfos[maxTasks];
	PrimInfo rinfos[maxTasks];
	
	/*! output data */
      private:
	List& lprims_o; 
	PrimInfo& linfo_o;
	List& rprims_o;
	PrimInfo& rinfo_o;
      };
    };
  }
}

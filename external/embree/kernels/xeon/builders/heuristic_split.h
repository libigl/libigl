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

#include "heuristic_object_partition_unaligned.h"
#include "heuristic_object_partition.h"
#include "heuristic_spatial_split.h"
#include "heuristic_strand_partition.h"
#include "heuristic_fallback.h"

namespace embree
{
  namespace isa
  {
    /*! Union type that selects between different splits */
    class __aligned(16) Split
    {
    private:
      
      enum Ty { OBJECT_SPLIT, OBJECT_SPLIT_UNALIGNED, SPATIAL_SPLIT, STRAND_SPLIT, FALLBACK_SPLIT };
      
      enum { 
	S0 = sizeof(ObjectPartitionUnaligned::Split), 
	S1 = sizeof(ObjectPartition::Split), 
	S2 = sizeof(SpatialSplit::Split), 
	S3 = sizeof(StrandSplit::Split),
	S01 = S0 > S1 ? S0 : S1,
	S23 = S2 > S3 ? S2 : S3,
	SIZE = S01 > S23 ? S01 : S23
      };
      
      typedef atomic_set<PrimRefBlockT<BezierPrim> > BezierRefList; //!< list of bezier primitives
      typedef atomic_set<PrimRefBlockT<PrimRef> > PrimRefList;   //!< list of general primitives
      
    public:
      
      /*! construct fallback split by default */
      __forceinline Split (bool isAligned = true) 
	: type(FALLBACK_SPLIT), sah(inf), isAligned(isAligned) {}
      
      /*! construction from object partitioning */
      __forceinline Split (ObjectPartitionUnaligned::Split& split) 
	: type(OBJECT_SPLIT_UNALIGNED), sah(split.splitSAH()), isAligned(false) { new (&data) ObjectPartitionUnaligned::Split(split); }
      
      /*! construction from object partitioning */
      __forceinline Split (const ObjectPartition::Split& split) 
	: type(OBJECT_SPLIT), sah(split.splitSAH()), isAligned(true) { new (&data) ObjectPartition::Split(split); }
      
      /*! construction from spatial split */
      __forceinline Split (SpatialSplit::Split& split) 
	: type(SPATIAL_SPLIT), sah(split.splitSAH()), isAligned(true) { new (&data) SpatialSplit::Split(split); }
      
      /*! construction from strand split */
      __forceinline Split (StrandSplit::Split& split) 
	: type(STRAND_SPLIT), sah(split.splitSAH()), isAligned(false) { new (&data) StrandSplit::Split(split); }
      
      /*! calculates surface area heuristic for performing the split */
      __forceinline float splitSAH() const { return sah; }
      
      /*! single threaded splitting into two sets */
      template<bool Parallel>
	void split(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimRefBlockAlloc<BezierPrim>& alloc, 
		   Scene* scene, BezierRefList& prims, 
		   BezierRefList& lprims_o, PrimInfo& linfo_o, 
		   BezierRefList& rprims_o, PrimInfo& rinfo_o) const
	{
	  switch (type) {
	  case OBJECT_SPLIT_UNALIGNED : ((ObjectPartitionUnaligned::Split*)&data)->split<Parallel>(threadIndex,threadCount,scheduler,alloc,prims,lprims_o,linfo_o,rprims_o,rinfo_o); break;
	  case OBJECT_SPLIT : ((ObjectPartition::Split*)&data)->split<Parallel>(threadIndex,threadCount,scheduler,alloc,prims,lprims_o,linfo_o,rprims_o,rinfo_o); break;
	  case SPATIAL_SPLIT: ((SpatialSplit::   Split*)&data)->split<Parallel>(threadIndex,threadCount,scheduler,alloc,scene,prims,lprims_o,linfo_o,rprims_o,rinfo_o); break;
	  case STRAND_SPLIT : ((StrandSplit::    Split*)&data)->split<Parallel>(threadIndex,threadCount,scheduler,alloc,prims,lprims_o,linfo_o,rprims_o,rinfo_o); break;
	  case FALLBACK_SPLIT: FallBackSplit::find(threadIndex,alloc,prims,lprims_o,linfo_o,rprims_o,rinfo_o); break;
	  default: THROW_RUNTIME_ERROR("internal error");
	  }
	}

      /*! single threaded splitting into two sets */
      template<bool Parallel>
	void split(size_t threadIndex, size_t threadCount, LockStepTaskScheduler* scheduler, PrimRefBlockAlloc<PrimRef>& alloc, 
		   Scene* scene, PrimRefList& prims, 
		   PrimRefList& lprims_o, PrimInfo& linfo_o, 
		   PrimRefList& rprims_o, PrimInfo& rinfo_o) const
	{
	  switch (type) {
	    //case OBJECT_SPLIT_UNALIGNED : ((ObjectPartitionUnaligned::Split*)&data)->split<Parallel>(threadIndex,threadCount,scheduler,alloc,prims,lprims_o,linfo_o,rprims_o,rinfo_o); break;
	  case OBJECT_SPLIT : ((ObjectPartition::Split*)&data)->split<Parallel>(threadIndex,threadCount,scheduler,alloc,prims,lprims_o,linfo_o,rprims_o,rinfo_o); break;
	  case SPATIAL_SPLIT: ((SpatialSplit::   Split*)&data)->split<Parallel>(threadIndex,threadCount,scheduler,alloc,scene,prims,lprims_o,linfo_o,rprims_o,rinfo_o); break;
	    //case STRAND_SPLIT : ((StrandSplit::    Split*)&data)->split<Parallel>(threadIndex,threadCount,scheduler,alloc,prims,lprims_o,linfo_o,rprims_o,rinfo_o); break;
	  case FALLBACK_SPLIT: FallBackSplit::find(threadIndex,alloc,prims,lprims_o,linfo_o,rprims_o,rinfo_o); break;
	  default: THROW_RUNTIME_ERROR("internal error");
	  }
	}
      
      private:
      __aligned(16) char data[SIZE]; //!< stores the different split types
      Ty type;                       //!< the type of split stored
      float sah;                     //!< the SAH of the stored split
    public:
      bool isAligned;
    };
  }
}

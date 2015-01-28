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

#ifndef __BVH4I_BUILDER_BINNING_H__
#define __BVH4I_BUILDER_BINNING_H__

#include "bvh4i.h"
#include "bvh4i_builder_util.h"
#include "bvh4i_builder_binner.h"

namespace embree
{
  namespace isa
  {

    template<int BINS>
    struct Mapping
    {
    public:
      
      __forceinline Mapping () {}
      
      __forceinline Mapping (const Centroid_Scene_AABB& bounds) 
        {
          const ssef centroidDiagonal = (ssef) bounds.centroid2.size();
          scale = select(centroidDiagonal != 0.0f,rcp(centroidDiagonal) * ssef(BINS * 0.99f),ssef(0.0f));
          ofs = (ssef) bounds.centroid2.lower;
        }
      
      /*! Computes the bin numbers for each dimension for a box. */
      __forceinline ssei bin_unsafe(const BBox3f& box) const {
        return floori((ssef(center2(box)) - ofs)*scale);
      }
      
      /*! Computes the bin numbers for each dimension for a box. */
      __forceinline ssei bin(const BBox3f& box) const {
#if defined (__SSE4_1__)
        return clamp(bin_unsafe(box),ssei(0),ssei(BINS-1));
#else
        ssei b = bin_unsafe(box);
        assert(b[0] >=0 && b[0] < BINS);
        assert(b[1] >=0 && b[1] < BINS);
        assert(b[2] >=0 && b[2] < BINS);
        return b;
#endif
      }
      
    public:
      ssef ofs;        //!< offset to compute bin
      ssef scale;      //!< scaling factor to compute bin
    };
    
    struct Split 
    {
      __forceinline Split () 
        : dim(-1), pos(-1), numLeft(-1), cost(pos_inf) {}
      
      /*! stream output */
      friend std::ostream& operator<<(std::ostream& cout, const Split& split) {
        return cout << "Split { " << 
          "dim = " << split.dim << 
          ", pos = " << split.pos << 
          ", numLeft = " << split.numLeft <<
          ", sah = " << split.cost << "}";
      }
      
    public:
      int dim;
      int pos;
      int numLeft;
      float cost;
    };
    
    template<int BINS>
      class Binner
    {
    public:
      
      /*! reset the binner */
      __forceinline void reset()
      {
        for (size_t i=0;i<BINS;i++) 
        {
          bounds[i][0] = empty;
          bounds[i][1] = empty;
          bounds[i][2] = empty;
          counts[i] = 0;
        }
      }
      
      /*! bin an array of primitives */
      void bin(const PrimRef* __restrict__ const prims, const size_t begin, const size_t end, const Mapping<BINS>& mapping);
      
      /*! bin an array of primitives and copy to destination array */
      void bin_copy(const PrimRef* __restrict__ const prims, const size_t begin, const size_t end, const Mapping<BINS>& mapping, PrimRef* __restrict__ const dst);
      
      /*! merge multiple binning infos into one */
      static void reduce(const Binner binners[], size_t num, Binner& binner_o);
      
      /*! calculate the best possible split */
      void best(Split& split, const Mapping<BINS>& mapping);
      
      /* inplace partitioning of a list of primitives */
      void partition(PrimRef*__restrict__ const prims,
                     const size_t begin,
                     const size_t end,
                     const Split& split,
                     const Mapping<BINS>& mapping,
                     BuildRecord& left,
                     BuildRecord& right);
      
    public:
      BBox3f bounds[BINS][4];
      ssei   counts[BINS];
    };
    
    bool split_fallback(PrimRef * __restrict__ const primref, BuildRecord& current, BuildRecord& leftChild, BuildRecord& rightChild);
    
    template<int BINS>
      class ParallelBinner
      {
      public:
        
        /*! parallel binbing of an array of primitives */
        void bin(BuildRecord& current, const PrimRef* src, PrimRef* dst, const size_t threadID, const size_t numThreads);
        
        /*! calculate the best possible split */
        void best(Split& split);
        
        /* parallel partitioning of a list of primitives */
        void partition(const PrimRef* src, PrimRef* dst, 
                       Split& split, 
                       BuildRecord &leftChild,
                       BuildRecord &rightChild,
                       const size_t threadID, const size_t numThreads);
        
      private:
        TASK_FUNCTION(ParallelBinner,parallelBinning);
        TASK_FUNCTION(ParallelBinner,parallelPartition);
        
      public:
        BuildRecord rec;
        Centroid_Scene_AABB left;
        Centroid_Scene_AABB right;
        Mapping<BINS> mapping;
        Split split;
        const PrimRef* src;
        PrimRef* dst;
        __aligned(64) AlignedAtomicCounter32 lCounter;
        __aligned(64) AlignedAtomicCounter32 rCounter;
        Binner<BINS> bin16;
        __aligned(64) Binner<BINS> global_bin16[MAX_MIC_THREADS];
      };
  }
}
#endif
  

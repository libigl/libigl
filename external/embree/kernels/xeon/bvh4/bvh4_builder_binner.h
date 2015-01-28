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

#ifndef __BVH4_BUILDER_BINNING_H__
#define __BVH4_BUILDER_BINNING_H__

#include "bvh4.h"
#include "../bvh4i/bvh4i_builder_util.h"

namespace embree
{
  namespace isa
  {
    struct BuildRef
    {
    public:
      __forceinline BuildRef () {}
      
      __forceinline BuildRef (const BBox3f& bounds, BVH4::NodeRef node) 
        : lower(bounds.lower), upper(bounds.upper), node(node)
      {
        if (node.isLeaf())
          lower.w = 0.0f;
        else
          lower.w = area(this->bounds());
      }
      
      __forceinline BBox3f bounds () const {
        return BBox3f(lower,upper);
      }
      
      friend bool operator< (const BuildRef& a, const BuildRef& b) {
        return a.lower.w < b.lower.w;
      }
      
    public:
      Vec3fa lower;
      Vec3fa upper;
      BVH4::NodeRef node;
    };

    template<int BINS>
    struct Mapping2
    {
    public:
      
      __forceinline Mapping2 () {}
      
      __forceinline Mapping2 (const Centroid_Scene_AABB& bounds) 
        {
          /* for toplevel builder we have to take geometry bounds here */
          const ssef geometryDiagonal = 2.0f * (ssef) bounds.geometry.size();
          scale = select(geometryDiagonal != 0.0f,rcp(geometryDiagonal) * ssef(BINS * 0.99f),ssef(0.0f));
          ofs = 2.0f * (ssef) bounds.geometry.lower;
        }
      
      /*! Computes the bin numbers for each dimension for a box. */
      __forceinline ssei bin_unsafe(const BBox3f& box) const {
        return floori((ssef(center2(box)) - ofs)*scale);
      }
      
      /*! Computes the bin numbers for each dimension for a box. */
      __forceinline ssei bin(const BBox3f& box) const {
#if defined(__SSE4_1__)
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
    
    struct Split2 
    {
      __forceinline Split2 () 
        : dim(-1), pos(-1), numLeft(-1), cost(pos_inf) {}
      
      /*! stream output */
      friend std::ostream& operator<<(std::ostream& cout, const Split2& split) {
        return cout << "Split2 { " << 
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
      class Binner2
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
      void bin(const BuildRef* __restrict__ const prims, const size_t begin, const size_t end, const Mapping2<BINS>& mapping);
      
      /*! bin an array of primitives and copy to destination array */
      void bin_copy(const BuildRef* __restrict__ const prims, const size_t begin, const size_t end, const Mapping2<BINS>& mapping, BuildRef* __restrict__ const dst);
      
      /*! merge multiple binning infos into one */
      static void reduce(const Binner2 binners[], size_t num, Binner2& binner_o);
      
      /*! calculate the best possible split */
      void best(Split2& split, const Mapping2<BINS>& mapping);
      
      /* inplace partitioning of a list of primitives */
      void partition(BuildRef*__restrict__ const prims,
                     const size_t begin,
                     const size_t end,
                     const Split2& split,
                     const Mapping2<BINS>& mapping,
                     BuildRecord& left,
                     BuildRecord& right);
      
    public:
      BBox3f bounds[BINS][4];
      ssei   counts[BINS];
    };
    
    bool split_fallback2(BuildRef * __restrict__ const primref, BuildRecord& current, BuildRecord& leftChild, BuildRecord& rightChild);
    
    template<int BINS>
      class ParallelBinner2
      {
      public:
        
        /*! parallel binbing of an array of primitives */
        void bin(BuildRecord& current, const BuildRef* src, BuildRef* dst, const size_t threadID, const size_t numThreads);
        
        /*! calculate the best possible split */
        void best(Split2& split);
        
        /* parallel partitioning of a list of primitives */
        void partition(const BuildRef* src, BuildRef* dst, 
                       Split2& split, 
                       BuildRecord &leftChild,
                       BuildRecord &rightChild,
                       const size_t threadID, const size_t numThreads);
        
      private:
        TASK_RUN_FUNCTION(ParallelBinner2,task_parallelBinning);
        TASK_RUN_FUNCTION(ParallelBinner2,task_parallelPartition);
        
      public:
        BuildRecord rec;
        Centroid_Scene_AABB left;
        Centroid_Scene_AABB right;
        Mapping2<BINS> mapping;
        Split2 split;
        const BuildRef* src;
        BuildRef* dst;
        __aligned(64) AlignedAtomicCounter32 lCounter;
        __aligned(64) AlignedAtomicCounter32 rCounter;
        Binner2<BINS> bin16;
        __aligned(64) Binner2<BINS> global_bin16[MAX_MIC_THREADS]; // FIXME: hardcoded number of threads
      };
  };
}

#endif


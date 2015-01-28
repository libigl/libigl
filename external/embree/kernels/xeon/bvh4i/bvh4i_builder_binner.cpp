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

#include "bvh4i_builder_binner.h"

#define QBVH_BUILDER_LEAF_ITEM_THRESHOLD 4

namespace embree
{
  namespace isa
  {
    template<int BINS>
    void Binner<BINS>::bin(const PrimRef * __restrict__ const prims,
                           const size_t begin,
                           const size_t end,
                           const Mapping<BINS>& mapping)
    {
      reset();
      if (end-begin == 0) return;
      
      /* unrolled binning loop */
      size_t i; 
      for (i=begin; i<end-1; i+=2)
      {
        /*! map even and odd primitive to bin */
        const BBox3f prim0 = prims[i+0].bounds(); const ssei bin0 = mapping.bin(prim0);
        const BBox3f prim1 = prims[i+1].bounds(); const ssei bin1 = mapping.bin(prim1);
        
        /*! increase bounds for bins for even primitive */
        const int b00 = bin0[0]; counts[b00][0]++; bounds[b00][0].extend(prim0);
        const int b01 = bin0[1]; counts[b01][1]++; bounds[b01][1].extend(prim0);
        const int b02 = bin0[2]; counts[b02][2]++; bounds[b02][2].extend(prim0);
        
        /*! increase bounds of bins for odd primitive */
        const int b10 = bin1[0]; counts[b10][0]++; bounds[b10][0].extend(prim1);
        const int b11 = bin1[1]; counts[b11][1]++; bounds[b11][1].extend(prim1);
        const int b12 = bin1[2]; counts[b12][2]++; bounds[b12][2].extend(prim1);
      }
      
      /*! for uneven number of primitives */
      if (i < end)
      {
        /*! map primitive to bin */
        const BBox3f prim0 = prims[i].bounds(); const ssei bin0 = mapping.bin(prim0);
        
        /*! increase bounds of bins */
        const int b00 = bin0[0]; counts[b00][0]++; bounds[b00][0].extend(prim0);
        const int b01 = bin0[1]; counts[b01][1]++; bounds[b01][1].extend(prim0);
        const int b02 = bin0[2]; counts[b02][2]++; bounds[b02][2].extend(prim0);
      }
    }
    
    template<int BINS>
    void Binner<BINS>::bin_copy(const PrimRef* __restrict__ const prims,
                                const size_t begin,
                                const size_t end,
                                const Mapping<BINS>& mapping,
                                PrimRef* __restrict__ const dest)
    {
      reset();
      if (end-begin == 0) return;
      
      /* unrolled binning loop */
      size_t i; 
      for (i=begin; i<end-1; i+=2)
      {
        /*! map even and odd primitive to bin */
        const BBox3f prim0 = prims[i+0].bounds(); const ssei bin0 = mapping.bin(prim0);
        const BBox3f prim1 = prims[i+1].bounds(); const ssei bin1 = mapping.bin(prim1);
        
        /*! increase bounds for bins for even primitive */
        const int b00 = bin0[0]; counts[b00][0]++; bounds[b00][0].extend(prim0);
        const int b01 = bin0[1]; counts[b01][1]++; bounds[b01][1].extend(prim0);
        const int b02 = bin0[2]; counts[b02][2]++; bounds[b02][2].extend(prim0);
        
        /*! increase bounds of bins for odd primitive */
        const int b10 = bin1[0]; counts[b10][0]++; bounds[b10][0].extend(prim1);
        const int b11 = bin1[1]; counts[b11][1]++; bounds[b11][1].extend(prim1);
        const int b12 = bin1[2]; counts[b12][2]++; bounds[b12][2].extend(prim1);
        
        /*! copy to destination */
        dest[i+0] = prims[i+0];
        dest[i+1] = prims[i+1];
      }
      
      /*! for uneven number of primitives */
      if (i < end)
      {
        /*! map primitive to bin */
        const BBox3f prim0 = prims[i].bounds(); const ssei bin0 = mapping.bin(prim0);
        
        /*! increase bounds of bins */
        const int b00 = bin0[0]; counts[b00][0]++; bounds[b00][0].extend(prim0);
        const int b01 = bin0[1]; counts[b01][1]++; bounds[b01][1].extend(prim0);
        const int b02 = bin0[2]; counts[b02][2]++; bounds[b02][2].extend(prim0);
        
        /*! copy to destination */
        dest[i+0] = prims[i+0];
      }
    }
    
    template<int BINS>
    void Binner<BINS>::reduce(const Binner binners[], size_t num, Binner& binner_o)
    {
      binner_o = binners[0];
      for (size_t tid=1; tid<num; tid++) 
      {
        const Binner& binner = binners[tid];
        for (size_t bin=0; bin<BINS; bin++) 
        {
          binner_o.bounds[bin][0].extend(binner.bounds[bin][0]);
          binner_o.bounds[bin][1].extend(binner.bounds[bin][1]);
          binner_o.bounds[bin][2].extend(binner.bounds[bin][2]);
          binner_o.counts[bin] += binner.counts[bin];
        }
      }
    }
    
    template<int BINS>
    void Binner<BINS>::best(Split& split, const Mapping<BINS>& mapping)
    {
      ssef rAreas[BINS];
      ssei rCounts[BINS];
      
      /* sweep from right to left and compute parallel prefix of merged bounds */
      ssei count = 0; BBox3f bx = empty; BBox3f by = empty; BBox3f bz = empty;
      for (size_t i=BINS-1; i>0; i--)
      {
        count += counts[i];
        rCounts[i] = count;
        bx.extend(bounds[i][0]); rAreas[i][0] = area(bx);
        by.extend(bounds[i][1]); rAreas[i][1] = area(by);
        bz.extend(bounds[i][2]); rAreas[i][2] = area(bz);
      }
      
      /* sweep from left to right and compute SAH */
      ssei ii = 1; ssef bestSAH = pos_inf; ssei bestPos = 0; ssei bestLeft = 0;
      count = 0; bx = empty; by = empty; bz = empty;
      for (size_t i=1; i<BINS; i++, ii+=1)
      {
        count += counts[i-1];
        bx.extend(bounds[i-1][0]); float Ax = area(bx);
        by.extend(bounds[i-1][1]); float Ay = area(by);
        bz.extend(bounds[i-1][2]); float Az = area(bz);
        const ssef lArea = ssef(Ax,Ay,Az,Az);
        const ssef rArea = rAreas[i];
        const ssei lCount = (count     +ssei(3)) >> 2;
        const ssei rCount = (rCounts[i]+ssei(3)) >> 2;
        const ssef sah = lArea*ssef(ssei_t(lCount)) + rArea*ssef(ssei_t(rCount));
        bestPos = select(sah < bestSAH,ii ,bestPos);
        bestLeft= select(sah < bestSAH,count,bestLeft);
        bestSAH = select(sah < bestSAH,sah,bestSAH);
      }
      
      /* find best dimension */
      for (size_t dim=0; dim<3; dim++) 
      {
        /* ignore zero sized dimensions */
        if (unlikely(mapping.scale[dim] == 0.0f)) 
          continue;
        
        /* test if this is a better dimension */
        if (bestSAH[dim] < split.cost && bestPos[dim] != 0) {
          split.dim = dim;
          split.pos = bestPos[dim];
          split.cost = bestSAH[dim];
          split.numLeft = bestLeft[dim];
        }
      }
    }
    
    template<typename PrimRef>
    __forceinline bool lt_split(const PrimRef *__restrict__ const aabb,
                                const unsigned int dim,
                                const float &c,
                                const float &s,
                                const int bestSplit) // FIXME: has to be singed int!!!!!!!!!
    {
      const ssef b_min(aabb->lower[dim]);
      const ssef b_max(aabb->upper[dim]);
      const ssef centroid_2 = b_min + b_max;
      const ssei binID = floori((centroid_2 - c)*s);
      return extract<0>(binID) < bestSplit;    
    }
    
    
    template<typename PrimRef>
    __forceinline bool ge_split(const PrimRef *__restrict__ const aabb,
                                const unsigned int dim,
                                const float &c,
                                const float &s,
                                const int bestSplit) // FIXME: has to be singed int!!!!!!!!!
    {
      const ssef b_min(aabb->lower[dim]);
      const ssef b_max(aabb->upper[dim]);
      const ssef centroid_2 = b_min + b_max;
      const ssei binID = floori((centroid_2 - c)*s);
      return extract<0>(binID) >= bestSplit;    
    }
    
    template<int BINS>
    void Binner<BINS>::partition(PrimRef *__restrict__ const prims,
                                 const size_t begin,
                                 const size_t end,
                                 const Split& split,
                                 const Mapping<BINS>& mapping,
                                 BuildRecord& left,
                                 BuildRecord& right)
    {
      Centroid_Scene_AABB local_left; local_left.reset();
      Centroid_Scene_AABB local_right; local_right.reset();
      
      assert(begin <= end);
      PrimRef* l = prims + begin;
      PrimRef* r = prims + end - 1;
      
      const float c = mapping.ofs[split.dim];
      const float s = mapping.scale[split.dim];
      const int bestSplitDim = split.dim;
      const int bestSplit = split.pos;
      
      ssef left_centroidMinAABB = (ssef) local_left.centroid2.lower;
      ssef left_centroidMaxAABB = (ssef) local_left.centroid2.upper;
      ssef left_sceneMinAABB    = (ssef) local_left.geometry.lower;
      ssef left_sceneMaxAABB    = (ssef) local_left.geometry.upper;
      
      ssef right_centroidMinAABB = (ssef) local_right.centroid2.lower;
      ssef right_centroidMaxAABB = (ssef) local_right.centroid2.upper;
      ssef right_sceneMinAABB    = (ssef) local_right.geometry.lower;
      ssef right_sceneMaxAABB    = (ssef) local_right.geometry.upper;
      
      while(1)
      {
        while (likely(l <= r && lt_split(l,bestSplitDim,c,s,bestSplit))) {
          const ssef b_min = load4f((float*)&l->lower);
          const ssef b_max = load4f((float*)&l->upper);
          const ssef centroid2 = b_min+b_max;
          left_centroidMinAABB = min(left_centroidMinAABB,centroid2);
          left_centroidMaxAABB = max(left_centroidMaxAABB,centroid2);
          left_sceneMinAABB    = min(left_sceneMinAABB,b_min);
          left_sceneMaxAABB    = max(left_sceneMaxAABB,b_max);
          ++l;
        }
        while (likely(l <= r && ge_split(r,bestSplitDim,c,s,bestSplit))) {
          const ssef b_min = load4f((float*)&r->lower);
          const ssef b_max = load4f((float*)&r->upper);
          const ssef centroid2 = b_min+b_max;
          right_centroidMinAABB = min(right_centroidMinAABB,centroid2);
          right_centroidMaxAABB = max(right_centroidMaxAABB,centroid2);
          right_sceneMinAABB    = min(right_sceneMinAABB,b_min);
          right_sceneMaxAABB    = max(right_sceneMaxAABB,b_max);
          --r;
        }
        if (r<l) break;
        
        const ssef r_min = load4f((float*)&l->lower);
        const ssef r_max = load4f((float*)&l->upper);
        const ssef r_centroid2 = r_min+r_max;
        right_centroidMinAABB = min(right_centroidMinAABB,r_centroid2);
        right_centroidMaxAABB = max(right_centroidMaxAABB,r_centroid2);
        right_sceneMinAABB    = min(right_sceneMinAABB,r_min);
        right_sceneMaxAABB    = max(right_sceneMaxAABB,r_max);
        const ssef l_min = load4f((float*)&r->lower);
        const ssef l_max = load4f((float*)&r->upper);
        const ssef l_centroid2 = l_min+l_max;
        left_centroidMinAABB = min(left_centroidMinAABB,l_centroid2);
        left_centroidMaxAABB = max(left_centroidMaxAABB,l_centroid2);
        left_sceneMinAABB    = min(left_sceneMinAABB,l_min);
        left_sceneMaxAABB    = max(left_sceneMaxAABB,l_max);
        store4f((float*)&l->lower,l_min);
        store4f((float*)&l->upper,l_max);
        store4f((float*)&r->lower,r_min);
        store4f((float*)&r->upper,r_max);
        l++; r--;
      }
      
      local_left.centroid2.lower = (Vec3fa) left_centroidMinAABB;
      local_left.centroid2.upper = (Vec3fa) left_centroidMaxAABB;
      local_left.geometry.lower = (Vec3fa) left_sceneMinAABB;
      local_left.geometry.upper = (Vec3fa) left_sceneMaxAABB;
      
      local_right.centroid2.lower = (Vec3fa) right_centroidMinAABB;
      local_right.centroid2.upper = (Vec3fa) right_centroidMaxAABB;
      local_right.geometry.lower = (Vec3fa) right_sceneMinAABB;
      local_right.geometry.upper = (Vec3fa) right_sceneMaxAABB;
      
      unsigned int center = l - prims;
      left.init(local_left,begin,center);
      right.init(local_right,center,end);
      
      assert(area(left.bounds.geometry) >= 0.0f);
      assert(area(left.bounds.centroid2) >= 0.0f);
      assert(area(right.bounds.geometry) >= 0.0f);
      assert(area(right.bounds.centroid2) >= 0.0f);
      
      assert( prims + begin <= l && l <= prims + end);
      assert( prims + begin <= r && r <= prims + end);
      
      assert(l <= prims + end);
      assert(center == begin+split.numLeft);
    }
    
    bool split_fallback(PrimRef * __restrict__ const primref, BuildRecord& current, BuildRecord& leftChild, BuildRecord& rightChild)
    {
      const unsigned int center = (current.begin + current.end)/2;
      
      Centroid_Scene_AABB left; left.reset();
      for (size_t i=current.begin; i<center; i++)
        left.extend(primref[i].bounds());
      leftChild.init(left,current.begin,center);
      
      Centroid_Scene_AABB right; right.reset();
      for (size_t i=center; i<current.end; i++)
        right.extend(primref[i].bounds());	
      rightChild.init(right,center,current.end);
      
      return true;
    }
    
    // =======================================================================================================
    // =======================================================================================================
    // =======================================================================================================
    
    template<int BINS>
    void ParallelBinner<BINS>::parallelBinning(const size_t threadID, const size_t numThreads)
    {
      Binner<16>& bin16 = this->global_bin16[threadID];
      const size_t startID = this->rec.begin + (threadID+0)*this->rec.items()/numThreads;
      const size_t endID   = this->rec.begin + (threadID+1)*this->rec.items()/numThreads;
      
      bin16.bin_copy(src,startID,endID,this->mapping,dst);
    }
    
    template<int BINS>
    void ParallelBinner<BINS>::bin(BuildRecord& current, const PrimRef* src, PrimRef* dst, const size_t threadID, const size_t numThreads) 
    {
      rec = current;
      mapping = Mapping<BINS>(current.bounds);
      left.reset();
      right.reset();
      this->src = src;
      this->dst = dst;
      LockStepTaskScheduler::dispatchTask(task_parallelBinning, this, threadID, numThreads );
      
      /* reduce binning information from all threads */
      Binner<BINS>::reduce(global_bin16,numThreads,bin16);
    }
    
    template<int BINS>
    void ParallelBinner<BINS>::best(Split& split) {
      bin16.best(split,mapping);
    }
    
    template<int BINS>
    void ParallelBinner<BINS>::parallelPartition(const size_t threadID, const size_t numThreads)
    {
      const size_t startID = this->rec.begin + (threadID+0)*this->rec.items()/numThreads;
      const size_t endID   = this->rec.begin + (threadID+1)*this->rec.items()/numThreads;
      
      /* load binning function */
      const unsigned int splitPos = this->split.pos;
      const unsigned int splitDim = this->split.dim;
      const float centroidBase = this->mapping.ofs[splitDim];
      const float centroidScale = this->mapping.scale[splitDim];
      
      /* compute items per thread that go to the 'left' and to the 'right' */
      int lnum[BINS];
      lnum[0] = this->global_bin16[threadID].counts[0][splitDim];
      for (size_t i=1; i<BINS; i++)
        lnum[i] = lnum[i-1] + this->global_bin16[threadID].counts[i][splitDim];
      
      const unsigned int localNumLeft = lnum[splitPos-1];
      const unsigned int localNumRight = (endID-startID) - localNumLeft;
      
      const unsigned int startLeft  = this->lCounter.add(localNumLeft);
      const unsigned int startRight = this->rCounter.add(localNumRight);
      
      PrimRef* __restrict__ src = (PrimRef*)this->src;
      PrimRef* __restrict__ dstLeft = dst + this->rec.begin + startLeft;
      PrimRef* __restrict__ dstRight = dst + this->rec.begin + startRight + this->split.numLeft;
      
      /* split into left and right */
      Centroid_Scene_AABB leftBounds; leftBounds.reset();
      Centroid_Scene_AABB rightBounds; rightBounds.reset();
      
      for (size_t i=startID; i<endID; i++)
      {
        if (likely(lt_split(&src[i],splitDim,centroidBase,centroidScale,splitPos))) {
          leftBounds.extend(src[i].bounds()); 
          *dstLeft++ = src[i];
        } else {
          rightBounds.extend(src[i].bounds()); 
          *dstRight++ = src[i];
        }
      }
      
      this->left .extend_atomic(leftBounds); 
      this->right.extend_atomic(rightBounds);  
    }
    
    template<int BINS>
    void ParallelBinner<BINS>::partition(const PrimRef* src, PrimRef* dst, 
                                         Split& split_i, 
                                         BuildRecord &leftChild,
                                         BuildRecord &rightChild,
                                         const size_t threadID, const size_t numThreads)
    {
      split = split_i;
      left.reset(); lCounter.reset(0);
      right.reset(); rCounter.reset(0); 
      this->src = src;
      this->dst = dst;
      LockStepTaskScheduler::dispatchTask( task_parallelPartition, this, threadID, numThreads );
      unsigned center = rec.begin + split.numLeft;
      assert(lCounter == split.numLeft);
      assert(rCounter == rec.items() - lCounter);
      leftChild.init(left,rec.begin,center);
      rightChild.init(right,center,rec.end);
    }
    
    /* explicit template instantiation */
    template class Binner<16>;
    template class ParallelBinner<16>;
  }
}


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

#include "common/primref.h"

namespace embree
{


  __aligned(64) static const int identity[16]         = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 };
  __aligned(64) static const int reverse_identity[16] = { 15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0 };

  __forceinline mic_f reverse(const mic_f &a) 
  {
    return _mm512_permutev_ps(load16i(reverse_identity),a);
  }
  
  __forceinline mic_f prefix_area_rl(const mic_f min_x,
				     const mic_f min_y,
				     const mic_f min_z,
				     const mic_f max_x,
				     const mic_f max_y,
				     const mic_f max_z)
  {
    const mic_f r_min_x = prefix_min(reverse(min_x));
    const mic_f r_min_y = prefix_min(reverse(min_y));
    const mic_f r_min_z = prefix_min(reverse(min_z));
    const mic_f r_max_x = prefix_max(reverse(max_x));
    const mic_f r_max_y = prefix_max(reverse(max_y));
    const mic_f r_max_z = prefix_max(reverse(max_z));
    
    const mic_f dx = r_max_x - r_min_x;
    const mic_f dy = r_max_y - r_min_y;
    const mic_f dz = r_max_z - r_min_z;
    
    const mic_f area_rl = (dx*dy+dx*dz+dy*dz) * mic_f(2.0f);
    return reverse(shl1_zero_extend(area_rl));
  }

  __forceinline mic_f prefix_area_lr(const mic_f min_x,
				     const mic_f min_y,
				     const mic_f min_z,
				     const mic_f max_x,
				     const mic_f max_y,
				     const mic_f max_z)
  {
    const mic_f r_min_x = prefix_min(min_x);
    const mic_f r_min_y = prefix_min(min_y);
    const mic_f r_min_z = prefix_min(min_z);
    const mic_f r_max_x = prefix_max(max_x);
    const mic_f r_max_y = prefix_max(max_y);
    const mic_f r_max_z = prefix_max(max_z);
    
    const mic_f dx = r_max_x - r_min_x;
    const mic_f dy = r_max_y - r_min_y;
    const mic_f dz = r_max_z - r_min_z;
  
    const mic_f area_lr = (dx*dy+dx*dz+dy*dz) * mic_f(2.0f);
    return area_lr;
  }


  __forceinline mic_i prefix_count(const mic_i c)
  {
    return prefix_sum(c);
  }

  static __forceinline mic_m lt_split(const mic_f &b_min,
				      const mic_f &b_max,
				      const mic_m &dim_mask,
				      const mic_f &c,
				      const mic_f &s,
				      const mic_f &bestSplit_f)
  {
    const mic_f centroid_2 = b_min + b_max;
    const mic_f binID = (centroid_2 - c)*s;
    return lt(dim_mask,binID,bestSplit_f);    
  }

  static __forceinline mic_m ge_split(const mic_f &b_min,
				      const mic_f &b_max,
				      const mic_m &dim_mask,
				      const mic_f &c,
				      const mic_f &s,
				      const mic_f &bestSplit_f)
  {
    const mic_f centroid_2 = b_min + b_max;
    const mic_f binID = (centroid_2 - c)*s;
    return ge(dim_mask,binID,bestSplit_f);    
  }


  template<class Primitive>
  __forceinline void fastbin(const Primitive * __restrict__ const aabb,
			     const unsigned int thread_start,
			     const unsigned int thread_end,
			     const mic_f &centroidBoundsMin_2,
			     const mic_f &scale,
			     mic_f lArea[3],
			     mic_f rArea[3],
			     mic_i lNum[3])
  {

    const mic_f init_min = mic_f::inf();
    const mic_f init_max = mic_f::minus_inf();
    const mic_i zero     = mic_i::zero();

    mic_f min_x0,min_x1,min_x2;
    mic_f min_y0,min_y1,min_y2;
    mic_f min_z0,min_z1,min_z2;
    mic_f max_x0,max_x1,max_x2;
    mic_f max_y0,max_y1,max_y2;
    mic_f max_z0,max_z1,max_z2;
    mic_i count0,count1,count2;

    min_x0 = init_min;
    min_x1 = init_min;
    min_x2 = init_min;
    min_y0 = init_min;
    min_y1 = init_min;
    min_y2 = init_min;
    min_z0 = init_min;
    min_z1 = init_min;
    min_z2 = init_min;

    max_x0 = init_max;
    max_x1 = init_max;
    max_x2 = init_max;
    max_y0 = init_max;
    max_y1 = init_max;
    max_y2 = init_max;
    max_z0 = init_max;
    max_z1 = init_max;
    max_z2 = init_max;

    count0 = zero;
    count1 = zero;
    count2 = zero;

    unsigned int start = thread_start;
    unsigned int end   = thread_end;

    if (unlikely(start % 2 != 0))
      {
	const mic2f bounds = aabb[start].getBounds();
	const mic_f b_min = bounds.x;
	const mic_f b_max = bounds.y;

	const mic_f centroid_2 = b_min + b_max; 
	const mic_i binID = convert_uint32((centroid_2 - centroidBoundsMin_2)*scale);
	// ------------------------------------------------------------------------      
	assert(0 <= binID[0] && binID[0] < 16);
	assert(0 <= binID[1] && binID[1] < 16);
	assert(0 <= binID[2] && binID[2] < 16);
	// ------------------------------------------------------------------------      
	const mic_i id = load16i(identity);
	const mic_m m_update_x = eq(id,swAAAA(binID));
	const mic_m m_update_y = eq(id,swBBBB(binID));
	const mic_m m_update_z = eq(id,swCCCC(binID));
	// ------------------------------------------------------------------------      
	min_x0 = mask_min(m_update_x,min_x0,min_x0,swAAAA(b_min));
	min_y0 = mask_min(m_update_x,min_y0,min_y0,swBBBB(b_min));
	min_z0 = mask_min(m_update_x,min_z0,min_z0,swCCCC(b_min));
	// ------------------------------------------------------------------------      
	max_x0 = mask_max(m_update_x,max_x0,max_x0,swAAAA(b_max));
	max_y0 = mask_max(m_update_x,max_y0,max_y0,swBBBB(b_max));
	max_z0 = mask_max(m_update_x,max_z0,max_z0,swCCCC(b_max));
	// ------------------------------------------------------------------------
	min_x1 = mask_min(m_update_y,min_x1,min_x1,swAAAA(b_min));
	min_y1 = mask_min(m_update_y,min_y1,min_y1,swBBBB(b_min));
	min_z1 = mask_min(m_update_y,min_z1,min_z1,swCCCC(b_min));      
	// ------------------------------------------------------------------------      
	max_x1 = mask_max(m_update_y,max_x1,max_x1,swAAAA(b_max));
	max_y1 = mask_max(m_update_y,max_y1,max_y1,swBBBB(b_max));
	max_z1 = mask_max(m_update_y,max_z1,max_z1,swCCCC(b_max));
	// ------------------------------------------------------------------------
	min_x2 = mask_min(m_update_z,min_x2,min_x2,swAAAA(b_min));
	min_y2 = mask_min(m_update_z,min_y2,min_y2,swBBBB(b_min));
	min_z2 = mask_min(m_update_z,min_z2,min_z2,swCCCC(b_min));
	// ------------------------------------------------------------------------      
	max_x2 = mask_max(m_update_z,max_x2,max_x2,swAAAA(b_max));
	max_y2 = mask_max(m_update_z,max_y2,max_y2,swBBBB(b_max));
	max_z2 = mask_max(m_update_z,max_z2,max_z2,swCCCC(b_max));
	// ------------------------------------------------------------------------
	count0 = mask_add(m_update_x,count0,count0,mic_i::one());
	count1 = mask_add(m_update_y,count1,count1,mic_i::one());
	count2 = mask_add(m_update_z,count2,count2,mic_i::one());      
	start++;	
      }
    assert(start % 2 == 0);

    if (unlikely(end % 2 != 0))
      {
	const mic2f bounds = aabb[end-1].getBounds();
	const mic_f b_min = bounds.x;
	const mic_f b_max = bounds.y;

	const mic_f centroid_2 = b_min + b_max; 
	const mic_i binID = convert_uint32((centroid_2 - centroidBoundsMin_2)*scale);
	// ------------------------------------------------------------------------      
	assert(0 <= binID[0] && binID[0] < 16);
	assert(0 <= binID[1] && binID[1] < 16);
	assert(0 <= binID[2] && binID[2] < 16);
	// ------------------------------------------------------------------------      
	const mic_i id = load16i(identity);
	const mic_m m_update_x = eq(id,swAAAA(binID));
	const mic_m m_update_y = eq(id,swBBBB(binID));
	const mic_m m_update_z = eq(id,swCCCC(binID));
	// ------------------------------------------------------------------------      
	min_x0 = mask_min(m_update_x,min_x0,min_x0,swAAAA(b_min));
	min_y0 = mask_min(m_update_x,min_y0,min_y0,swBBBB(b_min));
	min_z0 = mask_min(m_update_x,min_z0,min_z0,swCCCC(b_min));
	// ------------------------------------------------------------------------      
	max_x0 = mask_max(m_update_x,max_x0,max_x0,swAAAA(b_max));
	max_y0 = mask_max(m_update_x,max_y0,max_y0,swBBBB(b_max));
	max_z0 = mask_max(m_update_x,max_z0,max_z0,swCCCC(b_max));
	// ------------------------------------------------------------------------
	min_x1 = mask_min(m_update_y,min_x1,min_x1,swAAAA(b_min));
	min_y1 = mask_min(m_update_y,min_y1,min_y1,swBBBB(b_min));
	min_z1 = mask_min(m_update_y,min_z1,min_z1,swCCCC(b_min));      
	// ------------------------------------------------------------------------      
	max_x1 = mask_max(m_update_y,max_x1,max_x1,swAAAA(b_max));
	max_y1 = mask_max(m_update_y,max_y1,max_y1,swBBBB(b_max));
	max_z1 = mask_max(m_update_y,max_z1,max_z1,swCCCC(b_max));
	// ------------------------------------------------------------------------
	min_x2 = mask_min(m_update_z,min_x2,min_x2,swAAAA(b_min));
	min_y2 = mask_min(m_update_z,min_y2,min_y2,swBBBB(b_min));
	min_z2 = mask_min(m_update_z,min_z2,min_z2,swCCCC(b_min));
	// ------------------------------------------------------------------------      
	max_x2 = mask_max(m_update_z,max_x2,max_x2,swAAAA(b_max));
	max_y2 = mask_max(m_update_z,max_y2,max_y2,swBBBB(b_max));
	max_z2 = mask_max(m_update_z,max_z2,max_z2,swCCCC(b_max));
	// ------------------------------------------------------------------------
	count0 = mask_add(m_update_x,count0,count0,mic_i::one());
	count1 = mask_add(m_update_y,count1,count1,mic_i::one());
	count2 = mask_add(m_update_z,count2,count2,mic_i::one());      
	end--;	
      }
    assert(end % 2 == 0);

    const Primitive * __restrict__ aptr = aabb + start;

    prefetch<PFHINT_NT>(aptr);
    prefetch<PFHINT_L2>(aptr+2);
    prefetch<PFHINT_L2>(aptr+4);
    prefetch<PFHINT_L2>(aptr+6);
    prefetch<PFHINT_L2>(aptr+8);

    for (size_t j = start;j < end;j+=2,aptr+=2)
      {
	prefetch<PFHINT_L1>(aptr+2);
	prefetch<PFHINT_L2>(aptr+12);
	
#pragma unroll(2)
	for (size_t i=0;i<2;i++)
	  {

	    const mic2f bounds = aptr[i].getBounds();
	    const mic_f b_min  = bounds.x;
	    const mic_f b_max  = bounds.y;

	    const mic_f centroid_2 = b_min + b_max;
	    const mic_i binID = convert_uint32((centroid_2 - centroidBoundsMin_2)*scale);

	    assert(0 <= binID[0] && binID[0] < 16);
	    assert(0 <= binID[1] && binID[1] < 16);
	    assert(0 <= binID[2] && binID[2] < 16);

	    const mic_i id = load16i(identity);
	    const mic_m m_update_x = eq(id,swAAAA(binID));
	    const mic_m m_update_y = eq(id,swBBBB(binID));
	    const mic_m m_update_z = eq(id,swCCCC(binID));

	    min_x0 = mask_min(m_update_x,min_x0,min_x0,swAAAA(b_min));
	    min_y0 = mask_min(m_update_x,min_y0,min_y0,swBBBB(b_min));
	    min_z0 = mask_min(m_update_x,min_z0,min_z0,swCCCC(b_min));
	    // ------------------------------------------------------------------------      
	    max_x0 = mask_max(m_update_x,max_x0,max_x0,swAAAA(b_max));
	    max_y0 = mask_max(m_update_x,max_y0,max_y0,swBBBB(b_max));
	    max_z0 = mask_max(m_update_x,max_z0,max_z0,swCCCC(b_max));
	    // ------------------------------------------------------------------------
	    min_x1 = mask_min(m_update_y,min_x1,min_x1,swAAAA(b_min));
	    min_y1 = mask_min(m_update_y,min_y1,min_y1,swBBBB(b_min));
	    min_z1 = mask_min(m_update_y,min_z1,min_z1,swCCCC(b_min));      
	    // ------------------------------------------------------------------------      
	    max_x1 = mask_max(m_update_y,max_x1,max_x1,swAAAA(b_max));
	    max_y1 = mask_max(m_update_y,max_y1,max_y1,swBBBB(b_max));
	    max_z1 = mask_max(m_update_y,max_z1,max_z1,swCCCC(b_max));
	    // ------------------------------------------------------------------------
	    min_x2 = mask_min(m_update_z,min_x2,min_x2,swAAAA(b_min));
	    min_y2 = mask_min(m_update_z,min_y2,min_y2,swBBBB(b_min));
	    min_z2 = mask_min(m_update_z,min_z2,min_z2,swCCCC(b_min));
	    // ------------------------------------------------------------------------      
	    max_x2 = mask_max(m_update_z,max_x2,max_x2,swAAAA(b_max));
	    max_y2 = mask_max(m_update_z,max_y2,max_y2,swBBBB(b_max));
	    max_z2 = mask_max(m_update_z,max_z2,max_z2,swCCCC(b_max));
	    // ------------------------------------------------------------------------
	    count0 = mask_add(m_update_x,count0,count0,mic_i::one());
	    count1 = mask_add(m_update_y,count1,count1,mic_i::one());
	    count2 = mask_add(m_update_z,count2,count2,mic_i::one());      
	  }
      }

    prefetch<PFHINT_L1EX>(&rArea[0]);
    prefetch<PFHINT_L1EX>(&lArea[0]);
    prefetch<PFHINT_L1EX>(&lNum[0]);
    rArea[0] = prefix_area_rl(min_x0,min_y0,min_z0,max_x0,max_y0,max_z0);
    lArea[0] = prefix_area_lr(min_x0,min_y0,min_z0,max_x0,max_y0,max_z0);
    lNum[0]  = prefix_count(count0);

    prefetch<PFHINT_L1EX>(&rArea[1]);
    prefetch<PFHINT_L1EX>(&lArea[1]);
    prefetch<PFHINT_L1EX>(&lNum[1]);
    rArea[1] = prefix_area_rl(min_x1,min_y1,min_z1,max_x1,max_y1,max_z1);
    lArea[1] = prefix_area_lr(min_x1,min_y1,min_z1,max_x1,max_y1,max_z1);
    lNum[1]  = prefix_count(count1);

    prefetch<PFHINT_L1EX>(&rArea[2]);
    prefetch<PFHINT_L1EX>(&lArea[2]);
    prefetch<PFHINT_L1EX>(&lNum[2]);
    rArea[2] = prefix_area_rl(min_x2,min_y2,min_z2,max_x2,max_y2,max_z2);
    lArea[2] = prefix_area_lr(min_x2,min_y2,min_z2,max_x2,max_y2,max_z2);
    lNum[2]  = prefix_count(count2);
  }



  template<class Primitive>
  __forceinline void fastbin_xfm(const Primitive * __restrict__ const aabb,
				 const mic3f &cmat,
				 const unsigned int thread_start,
				 const unsigned int thread_end,
				 const mic_f &centroidBoundsMin_2,
				 const mic_f &scale,
				 mic_f lArea[3],
				 mic_f rArea[3],
				 mic_i lNum[3])
  {

    const mic_f init_min = mic_f::inf();
    const mic_f init_max = mic_f::minus_inf();
    const mic_i zero     = mic_i::zero();

    mic_f min_x0,min_x1,min_x2;
    mic_f min_y0,min_y1,min_y2;
    mic_f min_z0,min_z1,min_z2;
    mic_f max_x0,max_x1,max_x2;
    mic_f max_y0,max_y1,max_y2;
    mic_f max_z0,max_z1,max_z2;
    mic_i count0,count1,count2;

    min_x0 = init_min;
    min_x1 = init_min;
    min_x2 = init_min;
    min_y0 = init_min;
    min_y1 = init_min;
    min_y2 = init_min;
    min_z0 = init_min;
    min_z1 = init_min;
    min_z2 = init_min;

    max_x0 = init_max;
    max_x1 = init_max;
    max_x2 = init_max;
    max_y0 = init_max;
    max_y1 = init_max;
    max_y2 = init_max;
    max_z0 = init_max;
    max_z1 = init_max;
    max_z2 = init_max;

    count0 = zero;
    count1 = zero;
    count2 = zero;

    unsigned int start = thread_start;
    unsigned int end   = thread_end;

    const Primitive * __restrict__ aptr = aabb + start;

    prefetch<PFHINT_NT>(aptr);
    prefetch<PFHINT_L2>(aptr+2);
    prefetch<PFHINT_L2>(aptr+4);
    prefetch<PFHINT_L2>(aptr+6);
    prefetch<PFHINT_L2>(aptr+8);

    const mic_f c0 = cmat.x;
    const mic_f c1 = cmat.y;
    const mic_f c2 = cmat.z;

    for (size_t j = start;j < end;j++,aptr++)
      {
	prefetch<PFHINT_L1>(aptr+2);
	prefetch<PFHINT_L2>(aptr+12);
	
	const mic2f bounds = aptr->getBounds(c0,c1,c2);

	const mic_f b_min  = bounds.x;
	const mic_f b_max  = bounds.y;

	const mic_f centroid_2 = b_min + b_max;
	const mic_i binID_noclamp = convert_uint32((centroid_2 - centroidBoundsMin_2)*scale);
	const mic_i binID = min(max(binID_noclamp,mic_i::zero()),mic_i(15));
	assert(0 <= binID[0] && binID[0] < 16);
	assert(0 <= binID[1] && binID[1] < 16);
	assert(0 <= binID[2] && binID[2] < 16);

	const mic_i id = load16i(identity);
	const mic_m m_update_x = eq(id,swAAAA(binID));
	const mic_m m_update_y = eq(id,swBBBB(binID));
	const mic_m m_update_z = eq(id,swCCCC(binID));

	min_x0 = mask_min(m_update_x,min_x0,min_x0,swAAAA(b_min));
	min_y0 = mask_min(m_update_x,min_y0,min_y0,swBBBB(b_min));
	min_z0 = mask_min(m_update_x,min_z0,min_z0,swCCCC(b_min));
	// ------------------------------------------------------------------------      
	max_x0 = mask_max(m_update_x,max_x0,max_x0,swAAAA(b_max));
	max_y0 = mask_max(m_update_x,max_y0,max_y0,swBBBB(b_max));
	max_z0 = mask_max(m_update_x,max_z0,max_z0,swCCCC(b_max));
	// ------------------------------------------------------------------------
	min_x1 = mask_min(m_update_y,min_x1,min_x1,swAAAA(b_min));
	min_y1 = mask_min(m_update_y,min_y1,min_y1,swBBBB(b_min));
	min_z1 = mask_min(m_update_y,min_z1,min_z1,swCCCC(b_min));      
	// ------------------------------------------------------------------------      
	max_x1 = mask_max(m_update_y,max_x1,max_x1,swAAAA(b_max));
	max_y1 = mask_max(m_update_y,max_y1,max_y1,swBBBB(b_max));
	max_z1 = mask_max(m_update_y,max_z1,max_z1,swCCCC(b_max));
	// ------------------------------------------------------------------------
	min_x2 = mask_min(m_update_z,min_x2,min_x2,swAAAA(b_min));
	min_y2 = mask_min(m_update_z,min_y2,min_y2,swBBBB(b_min));
	min_z2 = mask_min(m_update_z,min_z2,min_z2,swCCCC(b_min));
	// ------------------------------------------------------------------------      
	max_x2 = mask_max(m_update_z,max_x2,max_x2,swAAAA(b_max));
	max_y2 = mask_max(m_update_z,max_y2,max_y2,swBBBB(b_max));
	max_z2 = mask_max(m_update_z,max_z2,max_z2,swCCCC(b_max));
	// ------------------------------------------------------------------------
	count0 = mask_add(m_update_x,count0,count0,mic_i::one());
	count1 = mask_add(m_update_y,count1,count1,mic_i::one());
	count2 = mask_add(m_update_z,count2,count2,mic_i::one());      
      }

    prefetch<PFHINT_L1EX>(&rArea[0]);
    prefetch<PFHINT_L1EX>(&lArea[0]);
    prefetch<PFHINT_L1EX>(&lNum[0]);
    rArea[0] = prefix_area_rl(min_x0,min_y0,min_z0,max_x0,max_y0,max_z0);
    lArea[0] = prefix_area_lr(min_x0,min_y0,min_z0,max_x0,max_y0,max_z0);
    lNum[0]  = prefix_count(count0);

    prefetch<PFHINT_L1EX>(&rArea[1]);
    prefetch<PFHINT_L1EX>(&lArea[1]);
    prefetch<PFHINT_L1EX>(&lNum[1]);
    rArea[1] = prefix_area_rl(min_x1,min_y1,min_z1,max_x1,max_y1,max_z1);
    lArea[1] = prefix_area_lr(min_x1,min_y1,min_z1,max_x1,max_y1,max_z1);
    lNum[1]  = prefix_count(count1);

    prefetch<PFHINT_L1EX>(&rArea[2]);
    prefetch<PFHINT_L1EX>(&lArea[2]);
    prefetch<PFHINT_L1EX>(&lNum[2]);
    rArea[2] = prefix_area_rl(min_x2,min_y2,min_z2,max_x2,max_y2,max_z2);
    lArea[2] = prefix_area_lr(min_x2,min_y2,min_z2,max_x2,max_y2,max_z2);
    lNum[2]  = prefix_count(count2);
  }


  template<unsigned int DISTANCE, class Primitive>
    __forceinline unsigned int partitionPrimitives_xfm(Primitive *__restrict__ aabb,
						       const mic3f &cmat,
						       const unsigned int begin,
						       const unsigned int end,
						       const unsigned int bestSplit,
						       const unsigned int bestSplitDim,
						       const mic_f &centroidBoundsMin_2,
						       const mic_f &scale,
						       Centroid_Scene_AABB & local_left,
						       Centroid_Scene_AABB & local_right)
    {
      assert(begin <= end);

      Primitive *__restrict__ l = aabb + begin;
      Primitive *__restrict__ r = aabb + end;

      const mic_f c = mic_f(centroidBoundsMin_2[bestSplitDim]);
      const mic_f s = mic_f(scale[bestSplitDim]);

      mic_f left_centroidMinAABB = broadcast4to16f(&local_left.centroid2.lower);
      mic_f left_centroidMaxAABB = broadcast4to16f(&local_left.centroid2.upper);
      mic_f left_sceneMinAABB    = broadcast4to16f(&local_left.geometry.lower);
      mic_f left_sceneMaxAABB    = broadcast4to16f(&local_left.geometry.upper);

      mic_f right_centroidMinAABB = broadcast4to16f(&local_right.centroid2.lower);
      mic_f right_centroidMaxAABB = broadcast4to16f(&local_right.centroid2.upper);
      mic_f right_sceneMinAABB    = broadcast4to16f(&local_right.geometry.lower);
      mic_f right_sceneMaxAABB    = broadcast4to16f(&local_right.geometry.upper);

      const mic_f bestSplit_f = mic_f(bestSplit);

      const mic_m dim_mask = mic_m::shift1[bestSplitDim];

      const mic_f c0 = cmat.x;
      const mic_f c1 = cmat.y;
      const mic_f c2 = cmat.z;

      while(1)
	{
	  while (likely(l < r)) 
	    {
	      
	      const mic2f bounds = l->getBounds(c0,c1,c2);
	      const mic_f b_min  = bounds.x;
	      const mic_f b_max  = bounds.y;

	      prefetch<PFHINT_L1EX>(l+2);	  
	      if (unlikely(ge_split(b_min,b_max,dim_mask,c,s,bestSplit_f))) break;
	      prefetch<PFHINT_L2EX>(l + DISTANCE + 4);	  
	      const mic_f centroid2 = b_min+b_max;
	      left_centroidMinAABB = min(left_centroidMinAABB,centroid2);
	      left_centroidMaxAABB = max(left_centroidMaxAABB,centroid2);
	      left_sceneMinAABB    = min(left_sceneMinAABB,b_min);
	      left_sceneMaxAABB    = max(left_sceneMaxAABB,b_max);
	      ++l;
	    }
	  while (likely(l < r)) 
	    {
	      const mic2f bounds = r->getBounds(c0,c1,c2);
	      const mic_f b_min  = bounds.x;
	      const mic_f b_max  = bounds.y;

	      prefetch<PFHINT_L1EX>(r-2);	  
	      if (unlikely(lt_split(b_min,b_max,dim_mask,c,s,bestSplit_f))) break;
	      prefetch<PFHINT_L2EX>(r - DISTANCE - 4);
	      const mic_f centroid2 = b_min+b_max;
	      right_centroidMinAABB = min(right_centroidMinAABB,centroid2);
	      right_centroidMaxAABB = max(right_centroidMaxAABB,centroid2);
	      right_sceneMinAABB    = min(right_sceneMinAABB,b_min);
	      right_sceneMaxAABB    = max(right_sceneMaxAABB,b_max);
	      --r;
	    }

	  if (unlikely(l == r)) {
	    const mic2f bounds = r->getBounds(c0,c1,c2);
	    const mic_f b_min  = bounds.x;
	    const mic_f b_max  = bounds.y;

	    if ( ge_split(b_min,b_max,dim_mask,c,s,bestSplit_f))
	      {
		const mic_f centroid2 = b_min+b_max;
		right_centroidMinAABB = min(right_centroidMinAABB,centroid2);
		right_centroidMaxAABB = max(right_centroidMaxAABB,centroid2);
		right_sceneMinAABB    = min(right_sceneMinAABB,b_min);
		right_sceneMaxAABB    = max(right_sceneMaxAABB,b_max);
	      }
	    else 
	      l++; 
	    break;
	  }

	  xchg(*l,*r);
	}


      store4f(&local_left.centroid2.lower,left_centroidMinAABB);
      store4f(&local_left.centroid2.upper,left_centroidMaxAABB);
      store4f(&local_left.geometry.lower,left_sceneMinAABB);
      store4f(&local_left.geometry.upper,left_sceneMaxAABB);

      store4f(&local_right.centroid2.lower,right_centroidMinAABB);
      store4f(&local_right.centroid2.upper,right_centroidMaxAABB);
      store4f(&local_right.geometry.lower,right_sceneMinAABB);
      store4f(&local_right.geometry.upper,right_sceneMaxAABB);

      assert( aabb + begin <= l && l <= aabb + end);
      assert( aabb + begin <= r && r <= aabb + end);

      return l - (aabb + begin);
    }


  class __aligned(64) Bin16
  {
  public:
    mic_f min_x[3];
    mic_f min_y[3];
    mic_f min_z[3];
    mic_f max_x[3];
    mic_f max_y[3];
    mic_f max_z[3];
    mic_i count[3];
    mic_i thread_count[3];

    Bin16() {}

    __forceinline void prefetchL2()
    {
      const unsigned int size = sizeof(Bin16);
#pragma unroll(8)
      for (size_t i=0;i<size;i+=64)
	prefetch<PFHINT_L2>(((char*)this) + i);
	
    }
    __forceinline void prefetchL2EX()
    {
      prefetch<PFHINT_L2EX>(&min_x[0]);
      prefetch<PFHINT_L2EX>(&min_x[1]);
      prefetch<PFHINT_L2EX>(&min_x[2]);

      prefetch<PFHINT_L2EX>(&min_y[0]);
      prefetch<PFHINT_L2EX>(&min_y[1]);
      prefetch<PFHINT_L2EX>(&min_y[2]);

      prefetch<PFHINT_L2EX>(&min_z[0]);
      prefetch<PFHINT_L2EX>(&min_z[1]);
      prefetch<PFHINT_L2EX>(&min_z[2]);

      prefetch<PFHINT_L2EX>(&max_x[0]);
      prefetch<PFHINT_L2EX>(&max_x[1]);
      prefetch<PFHINT_L2EX>(&max_x[2]);

      prefetch<PFHINT_L2EX>(&max_y[0]);
      prefetch<PFHINT_L2EX>(&max_y[1]);
      prefetch<PFHINT_L2EX>(&max_y[2]);

      prefetch<PFHINT_L2EX>(&max_z[0]);
      prefetch<PFHINT_L2EX>(&max_z[1]);
      prefetch<PFHINT_L2EX>(&max_z[2]);

      prefetch<PFHINT_L2EX>(&count[0]);
      prefetch<PFHINT_L2EX>(&count[1]);
      prefetch<PFHINT_L2EX>(&count[2]);

      prefetch<PFHINT_L2EX>(&thread_count[0]);
      prefetch<PFHINT_L2EX>(&thread_count[1]);
      prefetch<PFHINT_L2EX>(&thread_count[2]);
    }


    __forceinline void reset()
    {
      const mic_f init_min = mic_f::inf();
      const mic_f init_max = mic_f::minus_inf();
      const mic_i zero     = mic_i::zero();

      min_x[0] = init_min;
      min_x[1] = init_min;
      min_x[2] = init_min;

      min_y[0] = init_min;
      min_y[1] = init_min;
      min_y[2] = init_min;

      min_z[0] = init_min;
      min_z[1] = init_min;
      min_z[2] = init_min;

      max_x[0] = init_max;
      max_x[1] = init_max;
      max_x[2] = init_max;

      max_y[0] = init_max;
      max_y[1] = init_max;
      max_y[2] = init_max;

      max_z[0] = init_max;
      max_z[1] = init_max;
      max_z[2] = init_max;

      count[0] = zero;
      count[1] = zero;
      count[2] = zero;
    }


    __forceinline void merge(const Bin16& b)
    {
#pragma unroll(3)
      for (unsigned int i=0;i<3;i++)
	{
	  min_x[i] = min(min_x[i],b.min_x[i]);
	  min_y[i] = min(min_y[i],b.min_y[i]);
	  min_z[i] = min(min_z[i],b.min_z[i]);

	  max_x[i] = max(max_x[i],b.max_x[i]);
	  max_y[i] = max(max_y[i],b.max_y[i]);
	  max_z[i] = max(max_z[i],b.max_z[i]);

	  count[i] += b.count[i];
	}

    } 

  };

  __forceinline bool operator==(const Bin16 &a, const Bin16 &b) { 
    mic_m mask = 0xffff;
#pragma unroll(3)
    for (unsigned int i=0;i<3;i++)
      {
	mask &= eq(a.min_x[i],b.min_x[i]);
	mask &= eq(a.min_y[i],b.min_y[i]);
	mask &= eq(a.min_z[i],b.min_z[i]);

	mask &= eq(a.max_x[i],b.max_x[i]);
	mask &= eq(a.max_y[i],b.max_y[i]);
	mask &= eq(a.max_z[i],b.max_z[i]);

	mask &= eq(a.count[i],b.count[i]);
      }
    return mask == (mic_m)0xffff;
  };

  __forceinline bool operator!=(const Bin16 &a, const Bin16 &b) { 
    return !(a==b);
  }

  __forceinline std::ostream &operator<<(std::ostream &o, const Bin16 &v)
  {
#pragma unroll(3)
    for (unsigned int i=0;i<3;i++)
      {
	DBG_PRINT(v.min_x[i]);
	DBG_PRINT(v.min_y[i]);
	DBG_PRINT(v.min_z[i]);

	DBG_PRINT(v.max_x[i]);
	DBG_PRINT(v.max_y[i]);
	DBG_PRINT(v.max_z[i]);

	DBG_PRINT(v.count[i]);
      }

    return o;
  }



  template<class Primitive, bool NGO_OPTIMIZATION>
  __forceinline void fastbin_copy(const Primitive * __restrict__ const aabb,
				  Primitive * __restrict__ const tmp_aabb,
				  const unsigned int thread_start,
				  const unsigned int thread_end,
				  const mic_f &centroidBoundsMin_2,
				  const mic_f &scale,
				  Bin16 &bin16)
  {

    prefetch<PFHINT_NT>(&aabb[thread_start]);
    prefetch<PFHINT_NT>(&aabb[thread_end-1]);

    prefetch<PFHINT_L1EX>(&tmp_aabb[thread_start]);
    prefetch<PFHINT_L1EX>(&tmp_aabb[thread_end-1]);

    const mic_f init_min = mic_f::inf();
    const mic_f init_max = mic_f::minus_inf();
    const mic_i zero     = mic_i::zero();

    mic_f min_x0,min_x1,min_x2;
    mic_f min_y0,min_y1,min_y2;
    mic_f min_z0,min_z1,min_z2;
    mic_f max_x0,max_x1,max_x2;
    mic_f max_y0,max_y1,max_y2;
    mic_f max_z0,max_z1,max_z2;
    mic_i count0,count1,count2;

    min_x0 = init_min;
    min_x1 = init_min;
    min_x2 = init_min;
    min_y0 = init_min;
    min_y1 = init_min;
    min_y2 = init_min;
    min_z0 = init_min;
    min_z1 = init_min;
    min_z2 = init_min;

    max_x0 = init_max;
    max_x1 = init_max;
    max_x2 = init_max;
    max_y0 = init_max;
    max_y1 = init_max;
    max_y2 = init_max;
    max_z0 = init_max;
    max_z1 = init_max;
    max_z2 = init_max;

    count0 = zero;
    count1 = zero;
    count2 = zero;

    unsigned int start = thread_start;
    unsigned int end   = thread_end;

    if (start % 2 != 0)
      {
	const mic2f bounds = aabb[start].getBounds();
	const mic_f b_min  = bounds.x;
	const mic_f b_max  = bounds.y;

	const mic_f centroid_2 = b_min + b_max; 
	// const 
mic_i binID = convert_uint32((centroid_2 - centroidBoundsMin_2)*scale);
	// ------------------------------------------------------------------------      
        
        if (!(0 <= binID[0] && binID[0] < 16)) {
          std::cout << "binning assertion workaround" << std::endl;
          PRINT(binID);
          PRINT(bounds);
          PRINT(b_min);
          PRINT(b_max);
          PRINT(centroid_2);
          PRINT(centroidBoundsMin_2);
          PRINT(scale);
          binID = max(binID,mic_i(0));
          binID = min(binID,mic_i(15));
        } else {
          assert(0 <= binID[0] && binID[0] < 16);
          assert(0 <= binID[1] && binID[1] < 16);
          assert(0 <= binID[2] && binID[2] < 16);
        }
	// ------------------------------------------------------------------------      
	const mic_i id = load16i(identity);
	const mic_m m_update_x = eq(id,swAAAA(binID));
	const mic_m m_update_y = eq(id,swBBBB(binID));
	const mic_m m_update_z = eq(id,swCCCC(binID));
	// ------------------------------------------------------------------------      
	min_x0 = mask_min(m_update_x,min_x0,min_x0,swAAAA(b_min));
	min_y0 = mask_min(m_update_x,min_y0,min_y0,swBBBB(b_min));
	min_z0 = mask_min(m_update_x,min_z0,min_z0,swCCCC(b_min));
	// ------------------------------------------------------------------------      
	max_x0 = mask_max(m_update_x,max_x0,max_x0,swAAAA(b_max));
	max_y0 = mask_max(m_update_x,max_y0,max_y0,swBBBB(b_max));
	max_z0 = mask_max(m_update_x,max_z0,max_z0,swCCCC(b_max));
	// ------------------------------------------------------------------------
	min_x1 = mask_min(m_update_y,min_x1,min_x1,swAAAA(b_min));
	min_y1 = mask_min(m_update_y,min_y1,min_y1,swBBBB(b_min));
	min_z1 = mask_min(m_update_y,min_z1,min_z1,swCCCC(b_min));      
	// ------------------------------------------------------------------------      
	max_x1 = mask_max(m_update_y,max_x1,max_x1,swAAAA(b_max));
	max_y1 = mask_max(m_update_y,max_y1,max_y1,swBBBB(b_max));
	max_z1 = mask_max(m_update_y,max_z1,max_z1,swCCCC(b_max));
	// ------------------------------------------------------------------------
	min_x2 = mask_min(m_update_z,min_x2,min_x2,swAAAA(b_min));
	min_y2 = mask_min(m_update_z,min_y2,min_y2,swBBBB(b_min));
	min_z2 = mask_min(m_update_z,min_z2,min_z2,swCCCC(b_min));
	// ------------------------------------------------------------------------      
	max_x2 = mask_max(m_update_z,max_x2,max_x2,swAAAA(b_max));
	max_y2 = mask_max(m_update_z,max_y2,max_y2,swBBBB(b_max));
	max_z2 = mask_max(m_update_z,max_z2,max_z2,swCCCC(b_max));
	// ------------------------------------------------------------------------
	count0 = mask_add(m_update_x,count0,count0,mic_i::one());
	count1 = mask_add(m_update_y,count1,count1,mic_i::one());
	count2 = mask_add(m_update_z,count2,count2,mic_i::one());      
	tmp_aabb[start] = aabb[start];
	start++;	
      }
    assert(start % 2 == 0);

    if (end % 2 != 0)
      {
	const mic2f bounds = aabb[end-1].getBounds();
	const mic_f b_min  = bounds.x;
	const mic_f b_max  = bounds.y;

	const mic_f centroid_2 = b_min + b_max; 
	const mic_i binID = convert_uint32((centroid_2 - centroidBoundsMin_2)*scale);
	// ------------------------------------------------------------------------      
	assert(0 <= binID[0] && binID[0] < 16);
	assert(0 <= binID[1] && binID[1] < 16);
	assert(0 <= binID[2] && binID[2] < 16);
	// ------------------------------------------------------------------------      
	const mic_i id = load16i(identity);
	const mic_m m_update_x = eq(id,swAAAA(binID));
	const mic_m m_update_y = eq(id,swBBBB(binID));
	const mic_m m_update_z = eq(id,swCCCC(binID));
	// ------------------------------------------------------------------------      
	min_x0 = mask_min(m_update_x,min_x0,min_x0,swAAAA(b_min));
	min_y0 = mask_min(m_update_x,min_y0,min_y0,swBBBB(b_min));
	min_z0 = mask_min(m_update_x,min_z0,min_z0,swCCCC(b_min));
	// ------------------------------------------------------------------------      
	max_x0 = mask_max(m_update_x,max_x0,max_x0,swAAAA(b_max));
	max_y0 = mask_max(m_update_x,max_y0,max_y0,swBBBB(b_max));
	max_z0 = mask_max(m_update_x,max_z0,max_z0,swCCCC(b_max));
	// ------------------------------------------------------------------------
	min_x1 = mask_min(m_update_y,min_x1,min_x1,swAAAA(b_min));
	min_y1 = mask_min(m_update_y,min_y1,min_y1,swBBBB(b_min));
	min_z1 = mask_min(m_update_y,min_z1,min_z1,swCCCC(b_min));      
	// ------------------------------------------------------------------------      
	max_x1 = mask_max(m_update_y,max_x1,max_x1,swAAAA(b_max));
	max_y1 = mask_max(m_update_y,max_y1,max_y1,swBBBB(b_max));
	max_z1 = mask_max(m_update_y,max_z1,max_z1,swCCCC(b_max));
	// ------------------------------------------------------------------------
	min_x2 = mask_min(m_update_z,min_x2,min_x2,swAAAA(b_min));
	min_y2 = mask_min(m_update_z,min_y2,min_y2,swBBBB(b_min));
	min_z2 = mask_min(m_update_z,min_z2,min_z2,swCCCC(b_min));
	// ------------------------------------------------------------------------      
	max_x2 = mask_max(m_update_z,max_x2,max_x2,swAAAA(b_max));
	max_y2 = mask_max(m_update_z,max_y2,max_y2,swBBBB(b_max));
	max_z2 = mask_max(m_update_z,max_z2,max_z2,swCCCC(b_max));
	// ------------------------------------------------------------------------
	count0 = mask_add(m_update_x,count0,count0,mic_i::one());
	count1 = mask_add(m_update_y,count1,count1,mic_i::one());
	count2 = mask_add(m_update_z,count2,count2,mic_i::one());      
	tmp_aabb[end-1] = aabb[end-1];
	end--;	
      }
    assert(end % 2 == 0);

    const Primitive * __restrict__ aptr = aabb + start;

    prefetch<PFHINT_NT>(aptr);
    prefetch<PFHINT_L2>(aptr+2);
    prefetch<PFHINT_L2>(aptr+4);
    prefetch<PFHINT_L2>(aptr+6);
    prefetch<PFHINT_L2>(aptr+8);

    Primitive * __restrict__ tmp_aptr   = tmp_aabb + start;

    for (size_t j = start;j < end;j+=2,aptr+=2,tmp_aptr+=2)
      {
	const mic_f twoAABBs = load16f(aptr);

	prefetch<PFHINT_NT>(aptr+2);
	prefetch<PFHINT_L2>(aptr+12);

	if (NGO_OPTIMIZATION)
	  {
	    store16f_ngo(tmp_aptr,twoAABBs);
	  }
	else
	  {
	    tmp_aptr[0] = aptr[0]; 
	    tmp_aptr[1] = aptr[1]; 
	  }

#pragma unroll(2)
	for (size_t i=0;i<2;i++)
	  {
	    const mic2f bounds = aptr[i].getBounds();
	    const mic_f b_min  = bounds.x;
	    const mic_f b_max  = bounds.y;

	    const mic_f centroid_2 = b_min + b_max; 
	    const mic_i binID = convert_uint32((centroid_2 - centroidBoundsMin_2)*scale);

	    assert(0 <= binID[0] && binID[0] < 16); 
	    assert(0 <= binID[1] && binID[1] < 16); 
	    assert(0 <= binID[2] && binID[2] < 16); 

	    const mic_i id = load16i(identity);
	    const mic_m m_update_x = eq(id,swAAAA(binID));
	    const mic_m m_update_y = eq(id,swBBBB(binID));
	    const mic_m m_update_z = eq(id,swCCCC(binID));

	    min_x0 = mask_min(m_update_x,min_x0,min_x0,swAAAA(b_min));
	    min_y0 = mask_min(m_update_x,min_y0,min_y0,swBBBB(b_min));
	    min_z0 = mask_min(m_update_x,min_z0,min_z0,swCCCC(b_min));
	    // ------------------------------------------------------------------------      
	    max_x0 = mask_max(m_update_x,max_x0,max_x0,swAAAA(b_max));
	    max_y0 = mask_max(m_update_x,max_y0,max_y0,swBBBB(b_max));
	    max_z0 = mask_max(m_update_x,max_z0,max_z0,swCCCC(b_max));
	    // ------------------------------------------------------------------------
	    min_x1 = mask_min(m_update_y,min_x1,min_x1,swAAAA(b_min));
	    min_y1 = mask_min(m_update_y,min_y1,min_y1,swBBBB(b_min));
	    min_z1 = mask_min(m_update_y,min_z1,min_z1,swCCCC(b_min));      
	    // ------------------------------------------------------------------------      
	    max_x1 = mask_max(m_update_y,max_x1,max_x1,swAAAA(b_max));
	    max_y1 = mask_max(m_update_y,max_y1,max_y1,swBBBB(b_max));
	    max_z1 = mask_max(m_update_y,max_z1,max_z1,swCCCC(b_max));
	    // ------------------------------------------------------------------------
	    min_x2 = mask_min(m_update_z,min_x2,min_x2,swAAAA(b_min));
	    min_y2 = mask_min(m_update_z,min_y2,min_y2,swBBBB(b_min));
	    min_z2 = mask_min(m_update_z,min_z2,min_z2,swCCCC(b_min));
	    // ------------------------------------------------------------------------      
	    max_x2 = mask_max(m_update_z,max_x2,max_x2,swAAAA(b_max));
	    max_y2 = mask_max(m_update_z,max_y2,max_y2,swBBBB(b_max));
	    max_z2 = mask_max(m_update_z,max_z2,max_z2,swCCCC(b_max));
	    // ------------------------------------------------------------------------
	    count0 = mask_add(m_update_x,count0,count0,mic_i::one());
	    count1 = mask_add(m_update_y,count1,count1,mic_i::one());
	    count2 = mask_add(m_update_z,count2,count2,mic_i::one());      
	  }
      }

    bin16.prefetchL2EX();

    bin16.min_x[0] = min_x0;
    bin16.min_y[0] = min_y0;
    bin16.min_z[0] = min_z0;
    bin16.max_x[0] = max_x0;
    bin16.max_y[0] = max_y0;
    bin16.max_z[0] = max_z0;

    bin16.min_x[1] = min_x1;
    bin16.min_y[1] = min_y1;
    bin16.min_z[1] = min_z1;
    bin16.max_x[1] = max_x1;
    bin16.max_y[1] = max_y1;
    bin16.max_z[1] = max_z1;

    bin16.min_x[2] = min_x2;
    bin16.min_y[2] = min_y2;
    bin16.min_z[2] = min_z2;
    bin16.max_x[2] = max_x2;
    bin16.max_y[2] = max_y2;
    bin16.max_z[2] = max_z2;

    bin16.count[0] = count0;
    bin16.count[1] = count1;
    bin16.count[2] = count2;

    bin16.thread_count[0] = count0; 
    bin16.thread_count[1] = count1; 
    bin16.thread_count[2] = count2;     
  }

  template<class Primitive>
  static __forceinline mic_m lt_split(const Primitive *__restrict__ const aabb,
				      const unsigned int dim,
				      const mic_f &c,
				      const mic_f &s,
				      const mic_f &bestSplit_f)
  {
    const mic2f b = aabb->getBounds();
    const mic_f b_min = b.x;
    const mic_f b_max = b.y;
    const mic_f centroid_2 = b_min + b_max;
    const mic_f binID = (centroid_2 - c)*s;
    return lt(binID,bestSplit_f);    
  }

  template<class Primitive>
  static __forceinline mic_m ge_split(const Primitive *__restrict__ const aabb,
				      const unsigned int dim,
				      const mic_f &c,
				      const mic_f &s,
				      const mic_f &bestSplit_f)
  {
    const mic2f b = aabb->getBounds();
    const mic_f b_min = b.x;
    const mic_f b_max = b.y;
    const mic_f centroid_2 = b_min + b_max;
    const mic_f binID = (centroid_2 - c)*s;
    return ge(binID,bestSplit_f);    
  }



  template<unsigned int DISTANCE, class Primitive>
    __forceinline unsigned int partitionPrimitives(Primitive *__restrict__ aabb,
						   const unsigned int begin,
						   const unsigned int end,
						   const unsigned int bestSplit,
						   const unsigned int bestSplitDim,
						   const mic_f &centroidBoundsMin_2,
						   const mic_f &scale,
						   Centroid_Scene_AABB & local_left,
						   Centroid_Scene_AABB & local_right)
    {
      assert(begin <= end);

      Primitive *__restrict__ l = aabb + begin;
      Primitive *__restrict__ r = aabb + end;

      const mic_f c = mic_f(centroidBoundsMin_2[bestSplitDim]);
      const mic_f s = mic_f(scale[bestSplitDim]);

      mic_f left_centroidMinAABB = broadcast4to16f(&local_left.centroid2.lower);
      mic_f left_centroidMaxAABB = broadcast4to16f(&local_left.centroid2.upper);
      mic_f left_sceneMinAABB    = broadcast4to16f(&local_left.geometry.lower);
      mic_f left_sceneMaxAABB    = broadcast4to16f(&local_left.geometry.upper);

      mic_f right_centroidMinAABB = broadcast4to16f(&local_right.centroid2.lower);
      mic_f right_centroidMaxAABB = broadcast4to16f(&local_right.centroid2.upper);
      mic_f right_sceneMinAABB    = broadcast4to16f(&local_right.geometry.lower);
      mic_f right_sceneMaxAABB    = broadcast4to16f(&local_right.geometry.upper);

      const mic_f bestSplit_f = mic_f(bestSplit);

      const mic_m dim_mask = mic_m::shift1[bestSplitDim];

      while(1)
	{
	  while (likely(l < r)) 
	    {
	      
	      const mic2f bounds = l->getBounds();
	      const mic_f b_min  = bounds.x;
	      const mic_f b_max  = bounds.y;

	      prefetch<PFHINT_L1EX>(l+2);	  
	      if (unlikely(ge_split(b_min,b_max,dim_mask,c,s,bestSplit_f))) break;
	      prefetch<PFHINT_L2EX>(l + DISTANCE + 4);	  
	      const mic_f centroid2 = b_min+b_max;
	      left_centroidMinAABB = min(left_centroidMinAABB,centroid2);
	      left_centroidMaxAABB = max(left_centroidMaxAABB,centroid2);
	      left_sceneMinAABB    = min(left_sceneMinAABB,b_min);
	      left_sceneMaxAABB    = max(left_sceneMaxAABB,b_max);
	      ++l;
	    }
	  while (likely(l < r)) 
	    {
	      const mic2f bounds = r->getBounds();
	      const mic_f b_min  = bounds.x;
	      const mic_f b_max  = bounds.y;

	      prefetch<PFHINT_L1EX>(r-2);	  
	      if (unlikely(lt_split(b_min,b_max,dim_mask,c,s,bestSplit_f))) break;
	      prefetch<PFHINT_L2EX>(r - DISTANCE - 4);
	      const mic_f centroid2 = b_min+b_max;
	      right_centroidMinAABB = min(right_centroidMinAABB,centroid2);
	      right_centroidMaxAABB = max(right_centroidMaxAABB,centroid2);
	      right_sceneMinAABB    = min(right_sceneMinAABB,b_min);
	      right_sceneMaxAABB    = max(right_sceneMaxAABB,b_max);
	      --r;
	    }

	  if (unlikely(l == r)) {
	    const mic2f bounds = r->getBounds();
	    const mic_f b_min  = bounds.x;
	    const mic_f b_max  = bounds.y;

	    if ( ge_split(b_min,b_max,dim_mask,c,s,bestSplit_f))
	      {
		const mic_f centroid2 = b_min+b_max;
		right_centroidMinAABB = min(right_centroidMinAABB,centroid2);
		right_centroidMaxAABB = max(right_centroidMaxAABB,centroid2);
		right_sceneMinAABB    = min(right_sceneMinAABB,b_min);
		right_sceneMaxAABB    = max(right_sceneMaxAABB,b_max);
	      }
	    else 
	      l++; 
	    break;
	  }

	  xchg(*l,*r);
	}


      store4f(&local_left.centroid2.lower,left_centroidMinAABB);
      store4f(&local_left.centroid2.upper,left_centroidMaxAABB);
      store4f(&local_left.geometry.lower,left_sceneMinAABB);
      store4f(&local_left.geometry.upper,left_sceneMaxAABB);

      store4f(&local_right.centroid2.lower,right_centroidMinAABB);
      store4f(&local_right.centroid2.upper,right_centroidMaxAABB);
      store4f(&local_right.geometry.lower,right_sceneMinAABB);
      store4f(&local_right.geometry.upper,right_sceneMaxAABB);

      assert( aabb + begin <= l && l <= aabb + end);
      assert( aabb + begin <= r && r <= aabb + end);

      return l - (aabb + begin);
    }

  struct Split 
  {
    __forceinline void reset()
    {
      dim = -1;
      pos = -1;
      numLeft = -1;
      cost = pos_inf;
    }

    __forceinline Split () 
    {
      reset();
    }
    
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



    /* shared structure for multi-threaded binning and partitioning */
    struct __aligned(64) SharedBinningPartitionData
    {
      __aligned(64) BuildRecord rec;
      __aligned(64) Centroid_Scene_AABB left;
      __aligned(64) Centroid_Scene_AABB right;
      __aligned(64) Split split;
      __aligned(64) AlignedAtomicCounter32 lCounter;
      __aligned(64) AlignedAtomicCounter32 rCounter;
    };

};

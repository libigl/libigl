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

#ifndef __EMBREE_MIC_H__
#define __EMBREE_MIC_H__

#include "sys/platform.h"
#include "sys/intrinsics.h"

#include <immintrin.h>
#include <zmmintrin.h>

#define MIC_ALIGN __align(64)
#define _MM_SHUF_PERM(e3, e2, e1, e0) ((_MM_PERM_ENUM)((e3)*64 + (e2)*16 + (e1)*4 + (e0)))
#define _MM_SHUF_PERM_NONE _MM_SHUF_PERM(3,2,1,0)

namespace embree 
{
  class mic_m; 
  class mic_i; 
  class mic_f; 
}

#include "mic_m.h"
#include "mic_i.h"
#include "mic_f.h"

namespace embree 
{
  typedef mic_m micm_t;
  typedef mic_i mici_t;
  typedef mic_f micf_t;

  typedef mic_m micb_m;
  typedef mic_i mici_m;
  typedef mic_f micf_m;
}

#define PFHINT_L1   0
#define PFHINT_L2   1
#define PFHINT_NT   2
#define PFHINT_L1EX 3
#define PFHINT_L2EX 4
#define PFHINT_NTEX 5

namespace embree
{
  template<const unsigned int mode>
    __forceinline void prefetch(const void * __restrict__ const m)
  {
    if (mode == PFHINT_L1)
      _mm_prefetch((const char*)m,_MM_HINT_T0); 
    else if (mode == PFHINT_L2) 
      _mm_prefetch((const char*)m,_MM_HINT_T1); 
    else if (mode == PFHINT_NT) 
      _mm_prefetch((const char*)m,_MM_HINT_NTA); 
    else if (mode == PFHINT_L1EX)
      _mm_prefetch((const char*)m,_MM_HINT_ET0);  
    else if (mode == PFHINT_L2EX) 
      _mm_prefetch((const char*)m,_MM_HINT_ET2); 
    else if (mode == PFHINT_NTEX) 
      _mm_prefetch((const char*)m,_MM_HINT_ENTA); 
  }
  
  __forceinline void gather_prefetch(const mic_m &m_active,
                                     const void *const ptr, 			     
                                     const mic_i index, 
                                     const int mode = _MM_HINT_T2,
                                     const _MM_INDEX_SCALE_ENUM scale = _MM_SCALE_4,
                                     const _MM_UPCONV_PS_ENUM up = _MM_UPCONV_PS_NONE
                                     ) 
  {
    _mm512_mask_prefetch_i32extgather_ps(index,m_active,ptr,up,scale,mode);
  }
  
  __forceinline void scatter_prefetch(const mic_m &m_active,
                                      void *const ptr, 			     
                                      const mic_i index, 
                                      const int mode = _MM_HINT_ET2,
                                      const _MM_INDEX_SCALE_ENUM scale = _MM_SCALE_4,
                                      const _MM_UPCONV_PS_ENUM up = _MM_UPCONV_PS_NONE
                                      ) 
  {
    _mm512_mask_prefetch_i32extscatter_ps(ptr,m_active,index,up,scale,mode);
  }
  
  // ===========================
  // ====== prefetching ========
  // ===========================
  
  __forceinline  void evictL1(const void * __restrict__  m)
  {
    _mm_clevict(m,_MM_HINT_T0);
  }
  
  __forceinline  void evictL2(const void * __restrict__  m)
  {
    _mm_clevict(m,_MM_HINT_T1);
  }
  
  // ===============================
  // ====== dot, cross, ... ========
  // ===============================
  
  __forceinline mic_f lcross_zxy(const mic_f &ao, const mic_f &bo) {
    mic_f ao_bo = bo * swizzle(ao,_MM_SWIZ_REG_DACB);
    ao_bo = msub231(ao_bo,ao,swizzle(bo,_MM_SWIZ_REG_DACB));
    return ao_bo;
  }
  
  __forceinline mic_f ldot16_zxy(const mic_f &a,const mic_f &v0, const mic_f &v1, const mic_f &v2) 
  {
    mic_f v = v0 * swizzle(a,_MM_SWIZ_REG_BBBB);
    v = madd231(v,v1,swizzle(a,_MM_SWIZ_REG_CCCC));
    v = madd231(v,v2,swizzle(a,_MM_SWIZ_REG_AAAA));
    return v;
  }
  
  __forceinline mic_f ldot16_xyz(const mic_f &a,const mic_f &v0, const mic_f &v1, const mic_f &v2) 
  {
    mic_f v = v0 * swizzle(a,_MM_SWIZ_REG_AAAA);
    v = madd231(v, v1,swizzle(a,_MM_SWIZ_REG_BBBB));
    v = madd231(v, v2,swizzle(a,_MM_SWIZ_REG_CCCC));
    return v;
  }
  
  __forceinline mic_f lcross_xyz(const mic_f &a, const mic_f &b) 
  {
    mic_f c = b * swizzle(a,_MM_SWIZ_REG_DACB);
    c = msub231(c,a,swizzle(b,_MM_SWIZ_REG_DACB));
    c = swizzle(c,_MM_SWIZ_REG_DACB);
    return c;
  }
  
  __forceinline mic_f ldot3_xyz(const mic_f &ao, const mic_f &normal) 
  {
    mic_f vv = ao * normal;
    vv = _mm512_add_ps(vv,swizzle(vv,_MM_SWIZ_REG_CDAB));
    vv = _mm512_add_ps(vv,swizzle(vv,_MM_SWIZ_REG_BADC));
    return vv;        
  }

  __forceinline mic_f ldot3_xyz(const mic_m &m_mask, const mic_f &ao, const mic_f &normal) 
  {
    mic_f vv = _mm512_mask_mul_ps(ao,m_mask,ao,normal);
    vv = _mm512_add_ps(vv,swizzle(vv,_MM_SWIZ_REG_CDAB));
    vv = _mm512_add_ps(vv,swizzle(vv,_MM_SWIZ_REG_BADC));
    return vv;        
  }
  
  __forceinline mic_m lz_ldot3(const mic_f &ao, const mic_f &normal) 
  {
    mic_f vv = ao * normal;
    vv = _mm512_add_ps(vv,swizzle(vv,_MM_SWIZ_REG_CDAB));
    return lz_add(vv,swizzle(vv,_MM_SWIZ_REG_BADC));
  }
  
  __forceinline mic_m lz_ldot3(const mic_m &m_mask,const mic_f &ao, const mic_f &normal) 
  {
    mic_f vv = ao * normal;
    vv = _mm512_add_ps(vv,swizzle(vv,_MM_SWIZ_REG_CDAB));
    return lz_add(m_mask,vv,swizzle(vv,_MM_SWIZ_REG_BADC));
  }
  
  __forceinline mic_f ldot3_zxy(const mic_f &ao, const mic_f &normal) 
  {
    mic_f vv = ao * swizzle(normal,_MM_SWIZ_REG_DACB);
    vv = _mm512_add_ps(vv,swizzle(vv,_MM_SWIZ_REG_CDAB));
    vv = _mm512_add_ps(vv,swizzle(vv,_MM_SWIZ_REG_BADC));
    return vv;        
  }
  
  __forceinline mic_f lz_ldot3_zxy(const mic_f &ao, const mic_f &normal, mic_m &m_mask) 
  {
    mic_f vv = ao * swizzle(normal,_MM_SWIZ_REG_DACB);
    vv = _mm512_add_ps(vv,swizzle(vv,_MM_SWIZ_REG_CDAB));
    vv = lz_add(vv,swizzle(vv,_MM_SWIZ_REG_BADC),m_mask);
    return vv;        
  }
  
  __forceinline mic_f lz_ldot3(const mic_f &ao, const mic_f &normal, mic_m &m_mask) 
  {
    mic_f vv = ao * normal;
    vv = _mm512_add_ps(vv,swizzle(vv,_MM_SWIZ_REG_CDAB));
    vv = lz_add(vv,swizzle(vv,_MM_SWIZ_REG_BADC),m_mask);
    return vv;        
  }
  
  __forceinline mic_i cast_to_mic_i(mic_f f)
  {
    return _mm512_castps_si512(f);
  }
  
  __forceinline mic_f cast_to_mic_f(mic_i i)
  {
    return _mm512_castsi512_ps(i);
  }
  

  __forceinline mic_f reverse(const mic_f &a) 
  {
    return _mm512_permutev_ps(mic_i::reverse_identity(),a);
  }

  __forceinline mic_f SOAtoAOS_4f(const unsigned int index,
                                  const mic_f &x,
                                  const mic_f &y,
                                  const mic_f &z)
  {
    mic_f f = mic_f::zero();
    f = sel(0x1111,upconv1f((float*)&x + index),f);
    f = sel(0x2222,upconv1f((float*)&y + index),f);
    f = sel(0x4444,upconv1f((float*)&z + index),f);
    return f;
  }

  __forceinline mic_f gather_4f_zlc(const mic_i &v_mask,
                                    const void *__restrict__ const ptr0,
                                    const void *__restrict__ const ptr1,
                                    const void *__restrict__ const ptr2,
                                    const void *__restrict__ const ptr3) 
  {
    const mic_m m_00f0 = 0x00f0;
    const mic_m m_0f00 = 0x0f00;
    const mic_m m_f000 = 0xf000;
    
    mic_i v = v_mask & upconv4i((const int*)ptr0);
    v = mask_and(m_00f0,v,v_mask,upconv4i((const int*)ptr1));
    v = mask_and(m_0f00,v,v_mask,upconv4i((const int*)ptr2));
    v = mask_and(m_f000,v,v_mask,upconv4i((const int*)ptr3));
    return cast_to_mic_f(v);
  }

  __forceinline mic_f loadAOS(const float &x,const float &y, const float &z)
  {
    mic_f f = mic_f::zero();
    f = sel(0x1111,upconv1f(&x),f);
    f = sel(0x2222,upconv1f(&y),f);
    f = sel(0x4444,upconv1f(&z),f);
    return f;
  }
}

#endif

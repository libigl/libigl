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

#ifndef MIC_I_H
#define MIC_I_H

#include "sys/platform.h"

namespace embree
{
  /** mic_i class definition */
  class mic_i
  {
  private:
    __m512i v;
    
  public:
    
    /* construction */
    __forceinline mic_i() {};
    __forceinline mic_i(const int &i) { 
      v = _mm512_set_1to16_epi32(i);
      
    }
    __forceinline mic_i(const int &a, const int &b, const int &c, const int &d) { 
      v = _mm512_set_4to16_epi32(a,b,c,d);      
    }
    __forceinline mic_i(const mic_i &t) { v = t; };
    __forceinline mic_i(const __m512i &t) { v = t; };
    
    __forceinline explicit mic_i(const __m512 f) { 
      v = _mm512_cvtfxpnt_round_adjustps_epi32(f,_MM_FROUND_FLOOR,_MM_EXPADJ_NONE);
    }
    
    /* element access */  
    __forceinline int &operator[](const unsigned long i) { return ((int*)&v)[i]; };
    __forceinline const int &operator[](const unsigned long i) const { return ((int*)&v)[i]; };
    
    __forceinline operator __m512i () { return v; };
    __forceinline operator __m512i () const { return v; };
    
    __forceinline bool operator==(const mic_i &other) const { 
      return mic_m(_mm512_cmp_epi32_mask(v,other.v,_MM_CMPINT_EQ)) == MIC_M_ALL; 
    };
    
    /* assignment operators */
    __forceinline void operator=(const mic_i &f) { v = f.v; };
    __forceinline void operator=(const int &i) { 
      v = _mm512_set_1to16_epi32(i); 
    }
    
    __forceinline void operator+=(const mic_i &a) { 
      v = _mm512_add_epi32(v,a); 
    };
    
    __forceinline void operator-=(const mic_i &a) { 
      v = _mm512_sub_epi32(v,a);
    };
    
    __forceinline void operator*=(const mic_i &a) { 
      v = _mm512_mullo_epi32(v,a);
    };
    
    __forceinline void operator/=(const mic_i &a) { 
      v = _mm512_div_epi32(v,a); 
    }; 
    
    __forceinline void operator%=(const mic_i &a) { 
      v = _mm512_rem_epi32(v,a); 
    }; 
    
    __forceinline void operator&=(const mic_i &a) { 
      v = _mm512_and_epi32(v,a); 
    }; 
    
    __forceinline void operator|=(const mic_i &a) { 
      v = _mm512_or_epi32(v,a); 
    }; 
    
    __forceinline void operator^=(const mic_i &a) { 
      v = _mm512_xor_epi32(v,a); 
    }; 
    
    __forceinline void operator<<=(const mic_i &a) 
    { 
      v = _mm512_sllv_epi32(v,a); 
    }; 
    __forceinline void operator>>=(const mic_i &a) { 
      v = _mm512_srav_epi32(v,a); 
    }; 
    
    __forceinline void operator++() { 
      v = _mm512_add_epi32(v,mic_i::one()); 
    };
    
    __forceinline void operator--() { 
      v = _mm512_sub_epi32(v,mic_i::one()); 
    };
    
    /* constants */
    __forceinline static mic_i zero()      { 
      return _mm512_setzero_epi32(); 
    }
    __forceinline static mic_i one()       { return mic_i(1); }
    __forceinline static mic_i const1()    { return mic_i(1); }
    __forceinline static mic_i const2()    { return mic_i(2); }
    __forceinline static mic_i const3()    { return mic_i(3); }
    __forceinline static mic_i const4()    { return mic_i(4); }
    __forceinline static mic_i minus_one() { return mic_i(-1); }
    
    __forceinline static mic_i zlc4()      { 
      return _mm512_extload_epi32((void*)&int_zlc4,_MM_UPCONV_EPI32_NONE,_MM_BROADCAST_4X16,_MM_HINT_NONE);
    }
    
    static int int_identity[16];
    static int int_idMod4[16];
    static int int_idDiv4[16];
    static int int_pow4[16];
    static int int_zlc4[4];
    static int int_addtriID4[16];
    static unsigned int uint_shift1[32];
    static int int_aos2soa[16];
    static int int_reverse_identity[16];
    
    __forceinline static mic_i identity()  { 
      return _mm512_extload_epi32(int_identity,_MM_UPCONV_EPI32_NONE,_MM_BROADCAST_16X16,_MM_HINT_NONE);
    }; 
    
    __forceinline static mic_i idMod4()    { 
      return _mm512_extload_epi32(int_idMod4,_MM_UPCONV_EPI32_NONE,_MM_BROADCAST_16X16,_MM_HINT_NONE);
    }; 
    
    __forceinline static mic_i idDiv4()    { 
      return _mm512_extload_epi32(int_idDiv4,_MM_UPCONV_EPI32_NONE,_MM_BROADCAST_16X16,_MM_HINT_NONE);
    }; 
    __forceinline static mic_i addtriID4() { 
      return _mm512_extload_epi32(int_addtriID4,_MM_UPCONV_EPI32_NONE,_MM_BROADCAST_16X16,_MM_HINT_NONE);
    }; 
    __forceinline static mic_i aos2soa()   { 
      return _mm512_extload_epi32(int_aos2soa,_MM_UPCONV_EPI32_NONE,_MM_BROADCAST_16X16,_MM_HINT_NONE);
    }; 
    
    __forceinline static mic_i reverse_identity()  { 
      return _mm512_extload_epi32(int_reverse_identity,_MM_UPCONV_EPI32_NONE,_MM_BROADCAST_16X16,_MM_HINT_NONE);
    }; 
    
    __forceinline static unsigned int shift1(const unsigned int index) 
    { 
      assert(index < 32);
      return uint_shift1[index];
    }; 
    
  };
  
  __forceinline mic_i swizzle_i32(const mic_i &x, _MM_SWIZZLE_ENUM s) {
    return _mm512_swizzle_epi32(x,s);
  }
  
  __forceinline mic_i shuffle_i32(const mic_i &x,_MM_PERM_ENUM perm128, 
                                  _MM_SWIZZLE_ENUM perm32) 
  {
    return swizzle_i32(_mm512_permute4f128_epi32(x,perm128),perm32);
  }
  
  __forceinline mic_i shuffle_i32(const mic_m &mask, mic_i &v, const mic_i &x,_MM_PERM_ENUM perm128, 
                                  _MM_SWIZZLE_ENUM perm32) 
  {
    const __m512i p = _mm512_mask_permute4f128_epi32(v,mask,x,perm128);
    return _mm512_mask_swizzle_epi32(p,mask,x,perm32);  
  }
  
  __forceinline mic_i swAAAA_i(mic_i x) {
    return swizzle_i32(x,_MM_SWIZ_REG_AAAA);
  }
  
  __forceinline mic_i swBBBB_i(mic_i x) {
    return swizzle_i32(x,_MM_SWIZ_REG_BBBB);
  }
  
  __forceinline mic_i swCCCC_i(mic_i x) {
    return swizzle_i32(x,_MM_SWIZ_REG_CCCC);
  }
  
  __forceinline mic_i swDDDD_i(mic_i x) {
    return swizzle_i32(x,_MM_SWIZ_REG_DDDD);
  }
  
  /* loads */
  __forceinline mic_i upconv1i(const int *const ptr) { 
    return _mm512_extload_epi32(ptr,_MM_UPCONV_EPI32_NONE,_MM_BROADCAST_1X16,_MM_HINT_NONE);
  }
  
  __forceinline mic_i upconv1i_uint8(const unsigned char *const ptr) { 
    return _mm512_extload_epi32(ptr,_MM_UPCONV_EPI32_UINT8,_MM_BROADCAST_1X16,_MM_HINT_NONE);
  }
  
  __forceinline mic_i upconv4i(const int *const ptr) {
    return _mm512_extload_epi32(ptr,_MM_UPCONV_EPI32_NONE,_MM_BROADCAST_4X16,_MM_HINT_NONE);
  }
  
  __forceinline mic_i upconv16i(const int *const ptr) {
    return _mm512_extload_epi32(ptr,_MM_UPCONV_EPI32_NONE,_MM_BROADCAST_16X16,_MM_HINT_NONE);
  }
  
  __forceinline mic_i load16i(const int *const ptr) {
    return _mm512_extload_epi32(ptr,_MM_UPCONV_EPI32_NONE,_MM_BROADCAST_16X16,_MM_HINT_NONE);
    
  }
  
  __forceinline mic_i upconv16i_uint8(const unsigned char *const ptr) {
    return _mm512_extload_epi32(ptr,_MM_UPCONV_EPI32_UINT8,_MM_BROADCAST32_NONE,_MM_HINT_NONE);
  }
  
  __forceinline mic_i gather16i(const mic_m &mask, const int *const ptr, const mic_i &index,const _MM_INDEX_SCALE_ENUM scale) {
    return _mm512_mask_i32extgather_epi32(_mm512_undefined_epi32(),mask,index,ptr,_MM_UPCONV_EPI32_NONE,scale,0);
  }
  
  __forceinline mic_i gather16i(const mic_m &mask, mic_i &dest, const int *const ptr, const mic_i &index,const _MM_INDEX_SCALE_ENUM scale) {
    return _mm512_mask_i32extgather_epi32(dest,mask,index,ptr,_MM_UPCONV_EPI32_NONE,scale,0);
  }
  
  __forceinline void scatter16i(const mic_m &mask,int *const ptr, const mic_i &index,const mic_i &v, const _MM_INDEX_SCALE_ENUM scale) {
    _mm512_mask_i32extscatter_epi32((int*)ptr,mask,index,v,_MM_DOWNCONV_EPI32_NONE,scale,0);
  }
  
  /* stores */
  __forceinline void ustore16i(void *addr, const mic_i &reg) {
    _mm512_extpackstorelo_epi32((int*)addr+0  ,reg, _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
    _mm512_extpackstorehi_epi32((int*)addr+16 ,reg, _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
  }
  
  __forceinline void compactustore16i(const mic_m mask,void * addr, const mic_i &reg) {
    _mm512_mask_extpackstorelo_epi32((int*)addr+0  ,mask, reg, _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
    _mm512_mask_extpackstorehi_epi32((int*)addr+16 ,mask, reg, _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
  }
  
  __forceinline void compactustore16i_low(const mic_m mask, void *addr, const mic_i &reg) {
    _mm512_mask_extpackstorelo_epi32((int*)addr+0  ,mask, reg, _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
  }
  
  __forceinline void ustore16i_low(void *addr, const mic_i &reg) {
    _mm512_extpackstorelo_epi32((int*)addr+0  ,reg, _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
  }
  
  __forceinline void compactustore16i_high(const mic_m mask, int *addr, const mic_i &reg) {
    _mm512_mask_extpackstorehi_epi32((int*)addr+16  ,mask, reg, _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
  }
  
  // ------------------------------------------------------- 
  // comparisons 
  
  __forceinline mic_m eq(const mic_i &a, const mic_i &b)
  { 
    return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_EQ); 
  };
  
  __forceinline mic_m eq(const mic_m mask, const mic_i &a, const mic_i &b)
  { 
    return _mm512_mask_cmp_epi32_mask(mask,a,b,_MM_CMPINT_EQ); 
  };

  __forceinline mic_m lt(const mic_i &a, const mic_i &b)
  { 
    return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_LT); 
  };
  __forceinline mic_m lt(const mic_m mask, const mic_i &a, const mic_i &b)
  { 
    return _mm512_mask_cmp_epi32_mask(mask,a,b,_MM_CMPINT_LT); 
  };
  
  __forceinline mic_m le(const mic_i &a, const mic_i &b)
  { 
    return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_LE); 
  };
  
  __forceinline mic_m le(const mic_m mask, const mic_i &a, const mic_i &b)
  { 
    return _mm512_mask_cmp_epi32_mask(mask,a,b,_MM_CMPINT_LE); 
  };
  
  __forceinline mic_m ge(const mic_i &a, const mic_i &b)
  { 
    return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_GE); 
  };
  
  __forceinline mic_m ge(const mic_m mask, const mic_i &a, const mic_i &b)
  { 
    return _mm512_mask_cmp_epi32_mask(mask,a,b,_MM_CMPINT_GE); 
  };
  
  __forceinline mic_m gt(const mic_i &a, const mic_i &b)
  { 
    return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_GT); 
  };
  
  __forceinline mic_m gt(const mic_m mask, const mic_i &a, const mic_i &b)
  { 
    return _mm512_mask_cmp_epi32_mask(mask,a,b,_MM_CMPINT_GT); 
  };
  
  __forceinline mic_m ne(const mic_i &a, const mic_i &b)
  { 
    return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_NE); 
  };
  
  __forceinline mic_m ne(const mic_m mask, const mic_i &a, const mic_i &b)
  { 
    return _mm512_mask_cmp_epi32_mask(mask,a,b,_MM_CMPINT_NE); 
  };
  
  __forceinline mic_m operator!=(const mic_i &a, const mic_i &b) { 
    return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_NE); 
  };

  // ------------------------------------------------------- 
  // comparisons to zero
  
  __forceinline mic_m eqz(const mic_i &a) { return eq(a,mic_i::zero()); } 
  __forceinline mic_m eqz(const mic_m mask, const mic_i &a) { return eq(mask,a,mic_i::zero()); } 
  __forceinline mic_m lz(const mic_i &a) { return lt(a,mic_i::zero()); } 
  __forceinline mic_m lz(const mic_m mask, const mic_i &a) { return lt(mask,a,mic_i::zero()); } 
  __forceinline mic_m lez(const mic_i &a) { return le(a,mic_i::zero()); } 
  __forceinline mic_m lez(const mic_m mask, const mic_i &a) { return le(mask,a,mic_i::zero()); } 
  __forceinline mic_m gez(const mic_i &a) { return ge(a,mic_i::zero()); } 
  __forceinline mic_m gez(const mic_m mask, const mic_i &a) { return ge(mask,a,mic_i::zero()); } 
  __forceinline mic_m gz(const mic_i &a) { return gt(a,mic_i::zero()); } 
  __forceinline mic_m gz(const mic_m mask, const mic_i &a) { return gt(mask,a,mic_i::zero()); } 
  __forceinline mic_m nz(const mic_i &a) { return ne(a,mic_i::zero()); } 
  __forceinline mic_m nz(const mic_m mask, const mic_i &a) { return ne(mask,a,mic_i::zero()); } 
   
  // -------------------------------------------------------
  // arithmetic operators
  
  __forceinline mic_i operator+(const mic_i &a, const mic_i &b)
  { 
    return _mm512_add_epi32(a,b);
  };
  
  __forceinline mic_i operator-(const mic_i &a, const mic_i &b)
  { 
    return _mm512_sub_epi32(a,b);
  };
  
  __forceinline mic_i operator*(const mic_i &a, const mic_i &b)
  { 
    return _mm512_mullo_epi32(a,b);
  };
  
  __forceinline mic_i operator/(const mic_i &a, const mic_i &b)
  { 
    return _mm512_div_epi32(a,b);
  }; 
  
  __forceinline mic_i operator-(const mic_i &v)
  { 
    return mic_i::zero() - v; 
  };
  
  __forceinline mic_i operator+(const mic_i &v, const int &f)
  { 
    return v + mic_i(f); 
  };

  __forceinline mic_i operator+(const int &f,const mic_i &v)
  { 
    return v + mic_i(f); 
  };
  
  __forceinline mic_i operator-(const mic_i &v, const int &f)
  { 
    return v - mic_i(f);
  };
  
  __forceinline mic_i operator-(const int &f,const mic_i &v)
  { 
    return _mm512_subr_epi32(v,mic_i(f)); 
  };
  
  __forceinline mic_i operator*(const mic_i &v, const int &f) { 
    return v * mic_i(f); 
  };
  
  __forceinline mic_i operator*(const int &f,const mic_i &v) { 
    return v * mic_i(f);
  };
  
  __forceinline mic_i operator%(const mic_i &a, const mic_i &b)
  { 
    return _mm512_rem_epi32(a,b); 
  }; 
  
  __forceinline mic_i mask_add(const mic_m &mask, mic_i &c, const mic_i &a, const mic_i &b)
  { 
    return _mm512_mask_add_epi32(c,mask,a,b); 
  }; 
  
  __forceinline mic_i mask_sub(const mic_m &mask, mic_i &c, const mic_i &a, const mic_i &b)
  { 
    return _mm512_mask_sub_epi32(c,mask,a,b); 
  }; 
  
  __forceinline mic_i andn(const mic_i &a, const mic_i &b)
  { 
    return _mm512_andnot_epi32(a, b); 
  }
  
  // -------------------------------------------------------
  // bitwise operators
  
  __forceinline mic_i operator~(const mic_i &a) { return andn(a, mic_i(0xFFFFFFFF)); }
  
  __forceinline mic_i operator&(const mic_i &a, const mic_i &b) { 
  return _mm512_and_epi32(a,b);
  };
  
  __forceinline mic_i operator|(const mic_i &a, const mic_i &b) { 
    return _mm512_or_epi32(a,b);
  };
  
  __forceinline mic_i operator^(const mic_i &a, const mic_i &b) { 
    return _mm512_xor_epi32(a,b);
  };
  
  __forceinline mic_i operator<<(const mic_i &a, const mic_i &b) { 
    return _mm512_sllv_epi32(a,b);
  };
  
  __forceinline mic_i operator>>(const mic_i &a, const mic_i &b) { 
    return _mm512_srav_epi32(a,b);
  };
  
  /* reductions */
  __forceinline int reduce_add(mic_i a) { return _mm512_reduce_add_epi32(a); }
  __forceinline int reduce_mul(mic_i a) { return _mm512_reduce_mul_epi32(a); }
  __forceinline int reduce_min(mic_i a) { return _mm512_reduce_min_epi32(a); }
  __forceinline int reduce_max(mic_i a) { return _mm512_reduce_max_epi32(a); }
  __forceinline int reduce_and(mic_i a) { return _mm512_reduce_and_epi32(a); }
  
  // -------------------------------------------------------
  // selection
  
  __forceinline mic_i mask_and(const mic_m& m,mic_i &c, const mic_i &a, const mic_i &b)
  {
    return _mm512_mask_and_epi32(c,m,a,b);
  };
  
  __forceinline mic_i mask_or(const mic_m& m,mic_i &c, const mic_i &a, const mic_i &b)
  {
    return _mm512_mask_or_epi32(c,m,a,b);
  };
  
  __forceinline mic_i sel(const mic_m& m, const mic_i &a, const mic_i &b)
  { 
    return _mm512_mask_or_epi32(b,m,a,a); 
  };
  
  /* exchange */
  __forceinline void xchg(mic_m m, mic_i &a, mic_i &b) {
    mic_i tmp_a = a;
    a = sel(m,b,a);
    b = sel(m,tmp_a,b); 
  }
  
  /* min and max */
  __forceinline mic_i _min(mic_i a, mic_i b) { 
    return _mm512_min_epi32(a,b);
  }
  
  __forceinline mic_i _max(mic_i a, mic_i b) { 
    return _mm512_max_epi32(a,b);
  }
  
  
  /* horizontal min and max with broadcast */
  __forceinline mic_i set_min4(mic_i x) {
    x = _min(x,swizzle_i32(x,_MM_SWIZ_REG_BADC));
    x = _min(x,swizzle_i32(x,_MM_SWIZ_REG_CDAB));
    return x;
  }
  
  __forceinline mic_i set_max4(mic_i x) {
    x = _max(x,swizzle_i32(x,_MM_SWIZ_REG_BADC));
    x = _max(x,swizzle_i32(x,_MM_SWIZ_REG_CDAB));
    return x;
  }
  
  __forceinline mic_i set_min_lanes(mic_i x) {
    x = _min(x,_mm512_permute4f128_epi32(x, (_MM_PERM_ENUM)_MM_SHUF_PERM(2,3,0,1)));
    x = _min(x,_mm512_permute4f128_epi32(x, (_MM_PERM_ENUM)_MM_SHUF_PERM(1,0,3,2)));
    return x;
  }
  
  __forceinline mic_i set_max_lanes(mic_i x) {
    x = _max(x,_mm512_permute4f128_epi32(x, (_MM_PERM_ENUM)_MM_SHUF_PERM(2,3,0,1)));
    x = _max(x,_mm512_permute4f128_epi32(x, (_MM_PERM_ENUM)_MM_SHUF_PERM(1,0,3,2)));
    return x;
  }

  __forceinline mic_i set_min16(mic_i x) {
    return set_min_lanes(set_min4(x));
  }
  
  __forceinline mic_i set_max16(mic_i x) {
    return set_max_lanes(set_max4(x));
  }
  
  __forceinline mic_i gather16i_4i(const int *__restrict__ const ptr0,
                                   const int *__restrict__ const ptr1,
                                   const int *__restrict__ const ptr2,
                                   const int *__restrict__ const ptr3) 
  {
    mic_i v = upconv4i(ptr0);
    v = sel((mic_m)0xf0  ,upconv4i(ptr1),v);
    v = sel((mic_m)0xf00 ,upconv4i(ptr2),v);
    v = sel((mic_m)0xf000,upconv4i(ptr3),v);
    return v;
  }
  
  __forceinline void store16i(const mic_m &mask, void* __restrict__ addr, const mic_i& v2) {
    _mm512_mask_extstore_epi32(addr,mask,v2,_MM_DOWNCONV_EPI32_NONE,0);
  }
  
  __forceinline void store16i(void* __restrict__ addr, const mic_i& v2) {
    _mm512_extstore_epi32(addr,v2,_MM_DOWNCONV_EPI32_NONE,0);
  }
  
  __forceinline void store16i_uint8(const mic_m &mask, void* __restrict__ addr, const mic_i& v2) {
    _mm512_mask_extstore_epi32(addr,mask,v2,_MM_DOWNCONV_EPI32_UINT8,0);
  }
  
  __forceinline mic_i uload16i(const int *const addr) {
    mic_i r = _mm512_undefined_epi32();
    r =_mm512_extloadunpacklo_epi32(r, addr, _MM_UPCONV_EPI32_NONE, _MM_HINT_NONE);
    return _mm512_extloadunpackhi_epi32(r, addr+16, _MM_UPCONV_EPI32_NONE, _MM_HINT_NONE);  
  }
  
  __forceinline mic_i uload16i_low(const mic_m &mask, const void* addr) {
    mic_i v = _mm512_undefined_epi32();
    return _mm512_mask_extloadunpacklo_epi32(v, mask, addr, _MM_UPCONV_EPI32_NONE, _MM_HINT_NONE);
  }
  
  template<int D, int C, int B, int A>
    __forceinline mic_i lshuf(const mic_i &in)
  { 
    return _mm512_permute4f128_epi32(in,(_MM_PERM_ENUM)_MM_SHUF_PERM(D,C,B,A));
    
  }
  
  template<int D, int C, int B, int A>
    __forceinline mic_i lshuf(const mic_m &mask, mic_i &dest, const mic_i &in)
  { 
    return _mm512_mask_permute4f128_epi32(dest,mask,in,(_MM_PERM_ENUM)_MM_SHUF_PERM(D,C,B,A));
  }
  
  template<int lane>
    __forceinline mic_i shuffle_bc(const mic_i &x) {
    return shuffle_i32(x, (_MM_PERM_ENUM)_MM_SHUF_PERM(lane,lane,lane,lane),(_MM_SWIZZLE_ENUM)_MM_SHUF_PERM_NONE);
  }

  __forceinline mic_i prefix_sum(const mic_i &a)
  {
    mic_i v = a;
    v = mask_add(0xaaaa,v,v,swizzle_i32(v,_MM_SWIZ_REG_CDAB));
    v = mask_add(0xcccc,v,v,swizzle_i32(v,_MM_SWIZ_REG_BBBB));
    const mic_i shuf_v0 = shuffle_i32(v,(_MM_PERM_ENUM)_MM_SHUF_PERM(2,2,0,0),_MM_SWIZ_REG_DDDD);
    v = mask_add(0xf0f0,v,v,shuf_v0);
    const mic_i shuf_v1 = shuffle_i32(v,(_MM_PERM_ENUM)_MM_SHUF_PERM(1,1,0,0),_MM_SWIZ_REG_DDDD);
    v = mask_add(0xff00,v,v,shuf_v1);
    return v;  
  }
  
  __forceinline mic_i abs_diff(const mic_i &a, const mic_i &b)
  {
    return _max(a,b) - _min(a,b);
  }
  
  __forceinline mic_i set_or4(mic_i x) {
    x = x | swizzle_i32(x,_MM_SWIZ_REG_BADC);
    x = x | swizzle_i32(x,_MM_SWIZ_REG_CDAB);
    return x;
  }

  template<int i>
    __forceinline mic_i align_shift_right(const mic_i &a, const mic_i &b)
    {
      return _mm512_alignr_epi32(a,b,i); 
    };

  template<int i>
    __forceinline mic_i mask_align_shift_right(const mic_m &mask,mic_i &c,const mic_i &a, const mic_i &b)
    {
      return _mm512_mask_alignr_epi32(c,mask,a,b,i); 
    };
  
  __forceinline mic_m mask_test(const mic_m &m_mask, const mic_i &a, const mic_i &b) 
  {
    return _mm512_mask_test_epi32_mask(m_mask,a,b);
  }
  
  __forceinline void store16i_nr(void *__restrict__ ptr, const mic_i &a)
  {
    _mm512_storenr_ps(ptr,_mm512_castsi512_ps(a));
  }
  
  __forceinline void store16i_ngo(void *__restrict__ ptr, const mic_i &a)
  {
    _mm512_storenrngo_ps(ptr,_mm512_castsi512_ps(a));
  }
  
  __forceinline std::ostream &operator<<(std::ostream &o, const mic_i &v)
  {
    o << "[" << v[0];
    for (int i=1; i<16; i++) o << ", " << v[i];
    o << "]";
    return o;
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Select
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const mic_i select( const mic_m& s, const mic_i& t, const mic_i& f ) {
    return _mm512_mask_blend_epi32(s, f, t);
  }
}

#endif

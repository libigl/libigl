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

#ifndef MIC_F_H
#define MIC_F_H

#include "sys/platform.h"

namespace embree
{
  /** mic_f class definition */
  class mic_f 
  {
  public:
    __m512 v;
    
    /* construction */
    __forceinline mic_f() {};
    
    __forceinline mic_f(const float &f) { 
      v = _mm512_set_1to16_ps(f);
    }
    __forceinline mic_f(const float &a, const float &b, const float &c, const float &d) { 
      v = _mm512_set_4to16_ps(a,b,c,d);  
    }
    
    __forceinline mic_f(const __m512 &t) { v = t; };
    __forceinline mic_f(const mic_f &t) { v = t; };
    
    __forceinline explicit mic_f(const __m512i &a) { 
      v = _mm512_cvtfxpnt_round_adjustepi32_ps(a, _MM_FROUND_NO_EXC,_MM_EXPADJ_NONE);
    }
    __forceinline explicit mic_f(int *  a) { 
      __m512i i = _mm512_extload_epi32(a,_MM_UPCONV_EPI32_NONE,_MM_BROADCAST_1X16,_MM_HINT_NONE);
      v = _mm512_cvtfxpnt_round_adjustepi32_ps(i, _MM_FROUND_NO_EXC,_MM_EXPADJ_NONE);
    };
    
    __forceinline explicit mic_f(unsigned int *  a) { 
      __m512i i = _mm512_extload_epi32(a,_MM_UPCONV_EPI32_NONE,_MM_BROADCAST_1X16,_MM_HINT_NONE);
      v = _mm512_cvtfxpnt_round_adjustepi32_ps(i, _MM_FROUND_NO_EXC,_MM_EXPADJ_NONE);
    };
    
    /* element access */  
    __forceinline float &operator[](const unsigned long i) { return ((float*)&v)[i]; };
    __forceinline const float &operator[](const unsigned long i) const { return ((float*)&v)[i]; };
    
    __forceinline operator __m512 () { return v; };
    __forceinline operator __m512 () const { return v; };
    
    /* assignment operators */
    __forceinline void operator=(const mic_f &f) { v = f.v; };
    __forceinline void operator=(const float f) { v =  _mm512_set_1to16_ps(f);}
    
    __forceinline void operator+=(const mic_f &a) { v = _mm512_add_ps(v,a); };
    __forceinline void operator-=(const mic_f &a) { v = _mm512_sub_ps(v,a); };
    __forceinline void operator*=(const mic_f &a) { v = _mm512_mul_ps(v,a); };
    __forceinline void operator/=(const mic_f &a) { v = _mm512_div_ps(v,a); };
    
    __forceinline bool operator==(const mic_f &other) const { 
      // todo: test use instruction
      return mic_m(_mm512_cmp_ps_mask(v,other.v,_MM_CMPINT_EQ)) == MIC_M_ALL; 
    };
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline mic_f( ZeroTy   ) : v(_mm512_setzero()) {}
    __forceinline mic_f( OneTy    ) : v(_mm512_set_1to16_ps(1.0f)) {}
    __forceinline mic_f( PosInfTy ) : v(_mm512_set_1to16_ps(pos_inf)) {}
    __forceinline mic_f( NegInfTy ) : v(_mm512_set_1to16_ps(neg_inf)) {}
    __forceinline mic_f( StepTy   ) : v(_mm512_set_16to16_ps(15.0f,14.0f,13.0f,12.0f,11.0f,10.0f,9.0f,8.0f,7.0f,6.0f,5.0f,4.0f,3.0f,2.0f,1.0f,0.0f)) {}
    __forceinline mic_f( NaNTy    ) : v(_mm512_set_1to16_ps(nan)) {}
    
    /* constants */
    __forceinline static mic_f zero()   { return _mm512_setzero(); }
    __forceinline static mic_f one()    { return mic_f(1.0f); }
    __forceinline static mic_f const1() { return mic_f(1.0f); }
    __forceinline static mic_f const2() { return mic_f(2.0f); }
    __forceinline static mic_f const3() { return mic_f(3.0f); }
    __forceinline static mic_f const4() { return mic_f(4.0f); }
    
    __forceinline static mic_f minus_one() { return mic_f(-1.0f); }
    __forceinline static mic_f ulp() { return mic_f(float(embree::ulp)); }
    __forceinline static mic_f eps() { return mic_f(float(embree::ulp)); }
    __forceinline static mic_f inf() { return mic_f(float(pos_inf)); }
    __forceinline static mic_f minus_inf() { return mic_f(float(neg_inf)); }
    
    __forceinline static mic_f oneOver(const unsigned long index) { return mic_f(float_one_over[index]); }
    __forceinline static mic_f oneOver2()  { return mic_f(0.5f); }
    __forceinline static mic_f oneOver4()  { return mic_f(0.25f); }
    __forceinline static mic_f oneOver16() { return mic_f(0.0625f); }
    
    __forceinline static mic_f cancelLE() { 
      return _mm512_extload_ps(float_cancelLE,_MM_UPCONV_PS_NONE,_MM_BROADCAST_4X16,_MM_HINT_NONE);
    }
    
    static float float_idMod4[16];
    static float float_idDiv4[16];
    static float float_identity[16];
    static float float_cancelLE[4];
    static float float_one_over[32];
    
    __forceinline static mic_f idMod4() { return (mic_f&)float_idMod4[0]; }; 
    __forceinline static mic_f idDiv4() { return (mic_f&)float_idDiv4[0]; }; 
    __forceinline static mic_f identity() { return (mic_f&)float_identity[0]; }; 
    
  };

  
  __forceinline std::ostream &operator<<(std::ostream &o, const mic_f &v) 
  {
    o << "[" << v[0];
    for (int i=1; i<16; i++) o << ", " << v[i];
    o << "]";
    return o;
  } 
  
  __forceinline mic_f undefined() {
    return _mm512_undefined();
  }
  
  __forceinline mic_f load1f(const void *f) { 
    return _mm512_extload_ps(f,_MM_UPCONV_PS_NONE,_MM_BROADCAST_1X16,_MM_HINT_NONE);  
  }
  
  __forceinline mic_f load16f(const void *f) { 
    return _mm512_extload_ps(f,_MM_UPCONV_PS_NONE,_MM_BROADCAST_16X16,_MM_HINT_NONE);  
  }

  __forceinline mic_f load16f(const mic_f d, const mic_m valid, const void *f) { 
    return _mm512_mask_extload_ps(d,valid,f,_MM_UPCONV_PS_NONE,_MM_BROADCAST_16X16,_MM_HINT_NONE);  
  }
  
  __forceinline mic_f load4f(const void *f) { 
    return _mm512_extload_ps(f,_MM_UPCONV_PS_NONE,_MM_BROADCAST_4X16,0);  
  }
  
  __forceinline mic_f upconv1f(const float *const ptr) { 
    return _mm512_extload_ps(ptr,_MM_UPCONV_PS_NONE,_MM_BROADCAST_1X16,_MM_HINT_NONE);  
  }
  
  __forceinline mic_f upconv4f(const float *const ptr) {
    return _mm512_extload_ps(ptr,_MM_UPCONV_PS_NONE,_MM_BROADCAST_4X16,_MM_HINT_NONE);  
  }
  
  __forceinline mic_f upconv16f(const float *const ptr) {
    return _mm512_extload_ps(ptr,_MM_UPCONV_PS_NONE,_MM_BROADCAST_16X16,_MM_HINT_NONE);  
  }
  
  __forceinline mic_f upconv16_un8(const unsigned char *const ptr) {
    return _mm512_mul_ps(_mm512_extload_ps(ptr,_MM_UPCONV_PS_UINT8,_MM_BROADCAST_16X16,_MM_HINT_NONE),mic_f(1.0f/255.0f));  
  }
  
  /* loads */
  __forceinline mic_f uload16f(const float *const addr) {
    mic_f r = undefined();
    r =_mm512_extloadunpacklo_ps(r, addr, _MM_UPCONV_PS_NONE, _MM_HINT_NONE);
    return _mm512_extloadunpackhi_ps(r, addr+16, _MM_UPCONV_PS_NONE, _MM_HINT_NONE);  
  }
  
  __forceinline mic_f uload16f_low(const float *const addr) {
    mic_f r = undefined();
    return _mm512_extloadunpacklo_ps(r, addr, _MM_UPCONV_PS_NONE, _MM_HINT_NONE);
  }
  
  __forceinline mic_f uload16f_low(const mic_m &mask, const void* addr, const mic_f& v1) {
    return _mm512_mask_extloadunpacklo_ps(v1, mask, addr, _MM_UPCONV_PS_NONE, _MM_HINT_NONE);
  }
  
  __forceinline void ustore16f(float *addr, const mic_f &reg) {
    _mm512_extpackstorelo_ps(addr+0 ,reg, _MM_DOWNCONV_PS_NONE , 0);
    _mm512_extpackstorehi_ps(addr+16 ,reg, _MM_DOWNCONV_PS_NONE , 0);
  }
  
  __forceinline void compactustore16f(const mic_m &mask, float *addr, const mic_f &reg) {
    _mm512_mask_extpackstorelo_ps(addr+0 ,mask, reg, _MM_DOWNCONV_PS_NONE , 0);
    _mm512_mask_extpackstorehi_ps(addr+16 ,mask, reg, _MM_DOWNCONV_PS_NONE , 0);
  }
  
  __forceinline void compactustore16f_low(const mic_m &mask, float * addr, const mic_f &reg) {
    _mm512_mask_extpackstorelo_ps(addr+0 ,mask, reg, _MM_DOWNCONV_PS_NONE , 0);
  }
  
  __forceinline void ustore16f_low(float * addr, const mic_f &reg) {
    _mm512_extpackstorelo_ps(addr+0 ,reg, _MM_DOWNCONV_PS_NONE , 0);
  }
  
  __forceinline void compactustore16f_high(const mic_m &mask, float *addr, const mic_f &reg) {
    _mm512_extpackstorehi_ps(addr+0 ,reg, _MM_DOWNCONV_PS_NONE , 0);
  }
  
  __forceinline void store16f(const mic_m &mask, void* addr, const mic_f& v2) {
    _mm512_mask_extstore_ps(addr,mask,v2,_MM_DOWNCONV_PS_NONE,0);
  }
  
  __forceinline void store16f(void* addr, const mic_f& v2) {
    _mm512_extstore_ps(addr,v2,_MM_DOWNCONV_PS_NONE,0);
  }
  
  __forceinline void store4f(void* addr, const mic_f& v1) {
    assert((unsigned long)addr % 16 == 0);
    _mm512_mask_extpackstorelo_ps(addr,0xf, v1, _MM_DOWNCONV_PS_NONE , 0);
  }
  
  __forceinline void store16f_nr(void *__restrict__ ptr, const mic_f &a) {
    _mm512_storenr_ps(ptr,a);
  }
  
  __forceinline void store16f_ngo(void *__restrict__ ptr, const mic_f &a) {
    _mm512_storenrngo_ps(ptr,a);
  }
  
  __forceinline mic_f gather16f(const mic_m &mask, const float *const ptr, __m512i index, const _MM_INDEX_SCALE_ENUM scale) {
    mic_f r = undefined();
    return _mm512_mask_i32extgather_ps(r,mask,index,ptr,_MM_UPCONV_PS_NONE,scale,0);
  }
  
  __forceinline void scatter16f(const mic_m &mask,const float *const ptr, const __m512i index,const mic_f v, const _MM_INDEX_SCALE_ENUM scale) {
    _mm512_mask_i32extscatter_ps((void*)ptr,mask,index,v,_MM_DOWNCONV_PS_NONE,scale,0);
  }

  // ------------------------------------------------------- 
  // comparisons 
  
  __forceinline mic_m eq(const mic_f &a, const mic_f &b) { 
    return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_EQ); 
  }
  
  __forceinline mic_m eq(const mic_m &mask, const mic_f &a, const mic_f &b) { 
    return _mm512_mask_cmp_ps_mask(mask,a,b,_MM_CMPINT_EQ); 
  }
  
  __forceinline mic_m lt(const mic_f &a, const mic_f &b) { 
    return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_LT); 
  }
  
  __forceinline mic_m lt(const mic_m &mask, const mic_f &a, const mic_f &b) { 
    return _mm512_mask_cmp_ps_mask(mask,a,b,_MM_CMPINT_LT); 
  }
  
  __forceinline mic_m le(const mic_f &a, const mic_f &b) { 
    return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_LE); 
  }
  
  __forceinline mic_m le(const mic_m &mask, const mic_f &a, const mic_f &b) { 
    return _mm512_mask_cmp_ps_mask(mask,a,b,_MM_CMPINT_LE); 
  }
 
  __forceinline mic_m ge(const mic_f &a, const mic_f &b) { 
    return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_GE); 
  }

  __forceinline mic_m ge(const mic_m &mask, const mic_f &a, const mic_f &b) { 
    return _mm512_mask_cmp_ps_mask(mask,a,b,_MM_CMPINT_GE); 
  }
  
  __forceinline mic_m gt(const mic_f &a, const mic_f &b) { 
    return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_GT); 
  }

  __forceinline mic_m gt(const mic_m &mask, const mic_f &a, const mic_f &b) { 
    return _mm512_mask_cmp_ps_mask(mask,a,b,_MM_CMPINT_GT); 
  }
  
  __forceinline mic_m ne(const mic_f &a, const mic_f &b) { 
    return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_NE); 
  }
  
  __forceinline mic_m ne(const mic_m &mask, const mic_f &a, const mic_f &b) { 
    return _mm512_mask_cmp_ps_mask(mask,a,b,_MM_CMPINT_NE); 
  }

  __forceinline mic_m eqz(const mic_f &a) { return eq(a,mic_f::zero()); } 
  __forceinline mic_m eqz(const mic_m &mask, const mic_f &a) { return eq(mask,a,mic_f::zero()); } 
  __forceinline mic_m lz(const mic_f &a) { return lt(a,mic_f::zero()); } 
  __forceinline mic_m lz(const mic_m &mask, const mic_f &a) { return lt(mask,a,mic_f::zero()); } 
  __forceinline mic_m lez(const mic_f &a) { return le(a,mic_f::zero()); } 
  __forceinline mic_m lez(const mic_m &mask, const mic_f &a) { return le(mask,a,mic_f::zero()); } 
  __forceinline mic_m gez(const mic_f &a) { return ge(a,mic_f::zero()); } 
  __forceinline mic_m gez(const mic_m &mask, const mic_f &a) { return ge(mask,a,mic_f::zero()); } 
  __forceinline mic_m gz(const mic_f &a) { return gt(a,mic_f::zero()); } 
  __forceinline mic_m gz(const mic_m &mask, const mic_f &a) { return gt(mask,a,mic_f::zero()); } 
  __forceinline mic_m nz(const mic_f &a) { return ne(a,mic_f::zero()); } 
  __forceinline mic_m nz(const mic_m &mask, const mic_f &a) { return ne(mask,a,mic_f::zero()); } 

  // -------------------------------------------------------
  // arithmetic operators
  
  __forceinline mic_f operator+(const mic_f &a, const mic_f &b)
  { return _mm512_add_ps(a,b); };
  
  __forceinline mic_f operator-(const mic_f &a, const mic_f &b)
  { return _mm512_sub_ps(a,b); };
  
  __forceinline mic_f operator-(const mic_f &v) { 
    return _mm512_mul_ps(v,mic_f::minus_one()); 
  }
  
  __forceinline mic_f operator+(const mic_f &v, const float &f) { 
    return v + mic_f(f); 
  }
  
  __forceinline mic_f operator+(const float &f,const mic_f &v) { 
    return v + mic_f(f); 
  }
  
  __forceinline mic_f operator-(const mic_f &v, const float &f) { 
    return v - mic_f(f); 
  }
  
  __forceinline mic_f operator-(const float &f,const mic_f &v) { 
    return _mm512_subr_ps(v,mic_f(f)); 
  }
  
  __forceinline mic_f operator*(const mic_f &a, const mic_f &b) { 
    return _mm512_mul_ps(a,b); 
  }
  
  __forceinline mic_f operator*(const mic_f &v, const float &f) { 
    return v * mic_f(f); 
  }

  __forceinline mic_f operator*(const float &f,const mic_f &v) { 
    return v * mic_f(f); 
  }
  
  __forceinline mic_f operator/(const mic_f &a, const mic_f &b) 
  { return _mm512_div_ps(a,b); }; 
  
  __forceinline mic_f madd231(const mic_f &a, const mic_f &b, const mic_f &c) { 
    return _mm512_fmadd_ps(c,b,a);
  }
  
  __forceinline mic_f msub213(const mic_f &a, const mic_f &b, const mic_f &c) { 
    return _mm512_fmsub_ps(a,b,c);
  }
  
  __forceinline mic_f msub231(const mic_f &a, const mic_f &b, const mic_f &c) { 
    return _mm512_fmsub_ps(c,b,a);
  }
  
  __forceinline mic_f msubr231(const mic_f &a, const mic_f &b, const mic_f &c) { 
    return _mm512_fnmadd_ps(c,b,a);
  }

  __forceinline mic_f mul(const mic_m &mask,mic_f &c, const mic_f &a, const mic_f &b)
  { return _mm512_mask_mul_ps(c,mask,a,b); }; 
  
  __forceinline mic_f neg(const mic_m &mask,mic_f &c, const mic_f &a)
  { return _mm512_mask_addn_ps(c,mask,a,mic_f::zero()); }; 
  
  __forceinline mic_f mul(const mic_f &a, const mic_f &b)
  { return _mm512_mul_ps(a,b); }; 
  
  __forceinline mic_f add(const mic_f &a, const mic_f &b)
  { return _mm512_add_ps(a,b); }; 
  
  __forceinline mic_f subr(const mic_f &a,const mic_f &b)
  { return _mm512_subr_ps(b,a); };
  
  __forceinline mic_f subr_4f(const float *const a,const mic_f &b)
  { return _mm512_subr_ps(b,upconv4f(a)); };
  
  __forceinline mic_f subr_4f(const mic_m &mask,mic_f &old, const mic_f &a,const mic_f &b)
  { return _mm512_mask_subr_ps(old,mask,b,a); }
  
  __forceinline mic_f subr_4f(const mic_m &mask,mic_f &old,const float *const a,const mic_f &b)
  { return _mm512_mask_subr_ps(old,mask,b,upconv4f(a)); };
  
  __forceinline mic_f add_neg(const mic_f &a, const mic_f &b)
  {
    return  _mm512_add_ps(a,b) * mic_f::minus_one();  
  }
  
  __forceinline mic_f sel(const mic_m& m, const mic_f &a, const mic_f &b)
  {
    return _mm512_mask_mov_ps(b,m,a); 
  }; 
  
  __forceinline mic_f rcp(const mic_f &x)
  { 
    return _mm512_rcp23_ps(x);
  };
  
  __forceinline mic_f maxabs(const mic_f &a,const  mic_f &b) { return _mm512_maxabs_ps(a,b); }
  __forceinline mic_f abs(const mic_f &a) { 
    return _mm512_gmaxabs_ps(a,a); 
  }
  
  
  __forceinline mic_f sqr(const mic_f &a) { return a*a; } 
  __forceinline mic_f sqrt(const mic_f &a) { return _mm512_sqrt_ps(a); } 
  
  __forceinline mic_f length(const mic_f &a, const mic_f &b) { return _mm512_sqrt_ps(a*a+b*b); } 
  
  __forceinline mic_f rsqrt(const mic_f &a) {
    return _mm512_invsqrt_ps(a);
  }
  
  __forceinline mic_f vround(const mic_f &f, const int mode, const _MM_EXP_ADJ_ENUM exp = _MM_EXPADJ_NONE) { 
    return _mm512_round_ps(f,mode,exp); 
  }
  
  __forceinline mic_f fast_floor(const mic_f &a) { 
    return vround(a,(int)_MM_ROUND_MODE_DOWN);
  }
  
  __forceinline mic_f floor(const mic_f &a) { return _mm512_floor_ps(a); }
  __forceinline mic_f ceil(const mic_f &a) { return _mm512_ceil_ps(a); }
  __forceinline mic_f trunc(const mic_f &a) { return _mm512_trunc_ps(a); } 
  
  __forceinline mic_f exp(const mic_f &a) { return _mm512_exp_ps(a); }
  __forceinline mic_f exp2(const mic_f &a) { return _mm512_exp2_ps(a); }
  __forceinline mic_f pow(const mic_f &a, mic_f b) { return _mm512_pow_ps(a,b); }
  
  __forceinline mic_f log(const mic_f &a) { return _mm512_log_ps(a); }
  __forceinline mic_f log2(const mic_f &a) { return _mm512_log2_ps(a); }
  __forceinline mic_f log10(const mic_f &a) { return _mm512_log10_ps(a); }
  
  __forceinline mic_f sin(const mic_f &a) { return _mm512_sin_ps(a); } 
  __forceinline mic_f cos(const mic_f &a) { return _mm512_cos_ps(a); }
  __forceinline mic_f tan(const mic_f &a) { return _mm512_tan_ps(a); } 
  
  __forceinline mic_f asin(const mic_f &a) { return _mm512_asin_ps(a); }
  __forceinline mic_f acos(const mic_f &a) { return _mm512_acos_ps(a); }
  __forceinline mic_f atan(const mic_f &a) { return _mm512_atan_ps(a); }
  __forceinline mic_f atan2(const mic_f &a, mic_f b) { return _mm512_atan2_ps(a,b); }
  
  __forceinline mic_f sinh(const mic_f &a) { return _mm512_sinh_ps(a); } 
  __forceinline mic_f cosh(const mic_f &a) { return _mm512_cosh_ps(a); }
  __forceinline mic_f tanh(const mic_f &a) { return _mm512_tan_ps(a); } 
  
  __forceinline mic_f asinh(const mic_f &a) { return _mm512_asinh_ps(a); }
  __forceinline mic_f acosh(const mic_f &a) { return _mm512_acosh_ps(a); }
  __forceinline mic_f atanh(const mic_f &a) { return _mm512_atanh_ps(a); }
  
  __forceinline mic_f cbrt(const mic_f &a) { return _mm512_cbrt_ps(a); }
  __forceinline mic_f erf(const mic_f &a) { return _mm512_erf_ps(a); }
  __forceinline mic_f erfc(const mic_f &a) { return _mm512_erfc_ps(a); }
  __forceinline mic_f erfinv(const mic_f &a) { return _mm512_erfinv_ps(a); }
  __forceinline mic_f hypot(const mic_f &a, mic_f b) { return _mm512_hypot_ps(a,b); }
  __forceinline mic_f nearbyint(const mic_f &a) { return _mm512_nearbyint_ps(a); } 
  __forceinline mic_f rint(const mic_f &a) { return _mm512_rint_ps(a); } 
  
  /* reductions */
  __forceinline float reduce_add(const mic_f &a) { return _mm512_reduce_add_ps(a); }
  __forceinline float reduce_mul(const mic_f &a) { return _mm512_reduce_mul_ps(a); }
  __forceinline float reduce_min(const mic_f &a) { 
    return _mm512_reduce_min_ps(a); 
  }
  __forceinline float reduce_max(const mic_f &a) 
  { 
    return _mm512_reduce_max_ps(a); 
  }
  
  /* exchange */
  __forceinline void xchg(mic_m m, mic_f &a, mic_f &b) {
    mic_f tmp_a = a;
    a = sel(m,b,a);
    b = sel(m,tmp_a,b); 
  }
  
  __forceinline void xchg(mic_f &a, mic_f &b)
  {                                     
    const mic_f t = a;
    a = b;
    b = t;
  }                                                                     
  
  /* min and max */
  __forceinline mic_f _min(const mic_f &a, const mic_f &b) 
  { 
    return _mm512_gmin_ps(a,b); 
  }
  __forceinline mic_f _max(const mic_f &a, const mic_f &b) 
  { 
    return _mm512_gmax_ps(a,b); 
  }
  
  __forceinline mic_f _min(const mic_f &a, const mic_f &b,const mic_f &c, const mic_f &d) 
  {
    return _min(_min(a,b),_min(c,d));
  }
  
  __forceinline mic_f _max(const mic_f &a, const mic_f &b,const mic_f &c, const mic_f &d) 
  {
    return _max(_max(a,b),_max(c,d));
  }
  
  __forceinline mic_f swizzle(const mic_f &x, _MM_SWIZZLE_ENUM s) {
    return _mm512_swizzle_ps(x,s);
  }
  
  
  __forceinline mic_f shuffle(const mic_f &x, _MM_PERM_ENUM perm128, _MM_SWIZZLE_ENUM perm32) {
    return swizzle(_mm512_permute4f128_ps(x,perm128),perm32);
  }
  
  __forceinline mic_f shuffle(const mic_m &mask, mic_f &v, 
                              const mic_f &x,
                              _MM_PERM_ENUM perm128, 
                              _MM_SWIZZLE_ENUM perm32
                              ) {
    const __m512 p = _mm512_mask_permute4f128_ps(v,mask,x,perm128);
    return _mm512_mask_swizzle_ps(p,mask,x,perm32);  
  }
  
  template<int lane>
    __forceinline mic_f shuffle_bc(const mic_f &x) {
    return shuffle(x, _MM_SHUF_PERM(lane,lane,lane,lane), _MM_SHUF_PERM(3,2,1,0));
  }
  
  __forceinline mic_f swAAAA(const mic_f &x) {
  return swizzle(x,_MM_SWIZ_REG_AAAA);
  }
  
  __forceinline mic_f swBBBB(const mic_f &x) {
    return swizzle(x,_MM_SWIZ_REG_BBBB);
  }
  
  __forceinline mic_f swCCCC(const mic_f &x) {
    return swizzle(x,_MM_SWIZ_REG_CCCC);
  }
  
  __forceinline mic_f swDDDD(const mic_f &x) {
    return swizzle(x,_MM_SWIZ_REG_DDDD);
  }
  
  // =================================================
  // ===== horizontal min and max with broadcast ===== 
  // =================================================
  
  __forceinline mic_f set_min_lanes(mic_f x) {
    x = _min(x,_mm512_permute4f128_ps(x, (_MM_PERM_ENUM)_MM_SHUF_PERM(2,3,0,1)));
    x = _min(x,_mm512_permute4f128_ps(x, (_MM_PERM_ENUM)_MM_SHUF_PERM(1,0,3,2)));
    return x;
  }
  
  __forceinline mic_f set_max_lanes(mic_f x) {
    x = _max(x,_mm512_permute4f128_ps(x, (_MM_PERM_ENUM)_MM_SHUF_PERM(2,3,0,1)));
    x = _max(x,_mm512_permute4f128_ps(x, (_MM_PERM_ENUM)_MM_SHUF_PERM(1,0,3,2)));
    return x;
  }
  
  __forceinline mic_f set_min4(mic_f x) {
    x = _min(x,swizzle(x,_MM_SWIZ_REG_BADC));
    x = _min(x,swizzle(x,_MM_SWIZ_REG_CDAB));
    return x;
  }
  
  __forceinline mic_f set_max4(mic_f x) {
    x = _max(x,swizzle(x,_MM_SWIZ_REG_BADC));
    x = _max(x,swizzle(x,_MM_SWIZ_REG_CDAB));
    return x;
  }
  
  __forceinline mic_f set_add4(mic_f x) {
    x = _mm512_add_ps(x,swizzle(x,_MM_SWIZ_REG_BADC));
    x = _mm512_add_ps(x,swizzle(x,_MM_SWIZ_REG_CDAB));
    return x;
  }
  
  __forceinline mic_f set_min16(mic_f x) {
    return set_min_lanes(set_min4(x));
  }
  
  __forceinline mic_f set_max16(mic_f x) {
    return set_max_lanes(set_max4(x));
  }
  
  __forceinline mic_f gather16f_4f(const float *__restrict__ const ptr0,
                                   const float *__restrict__ const ptr1,
                                   const float *__restrict__ const ptr2,
                                   const float *__restrict__ const ptr3) 
  {
    mic_f v = upconv4f(ptr0);
    v = sel((mic_m)0x00f0,upconv4f(ptr1),v);
    v = sel((mic_m)0x0f00,upconv4f(ptr2),v);
    v = sel((mic_m)0xf000,upconv4f(ptr3),v);
    return v;
  }
  
  __forceinline mic_f gather16f_4f(const mic_f &v0,
                                   const mic_f &v1,
                                   const mic_f &v2,
                                   const mic_f &v3) 
  {
    mic_f v = v0;
    v = sel((mic_m)0xf0  ,v1,v);
    v = sel((mic_m)0xf00 ,v2,v);
    v = sel((mic_m)0xf000,v3,v);
    return v;
  }
  
  
  template<const int D, const int C, const int B, const int A> 
    __forceinline mic_f lshuf(const mic_f &in)
  { 
    return _mm512_permute4f128_ps(in,(_MM_PERM_ENUM)_MM_SHUF_PERM(D,C,B,A));
  }
  
  
  template<const int D, const int C, const int B, const int A> 
    __forceinline mic_f lshuf(const mic_m &mask, mic_f &dest, const mic_f &in)
  { 
    return _mm512_mask_permute4f128_ps(dest,mask,in,(_MM_PERM_ENUM) _MM_SHUF_PERM(D,C,B,A));
  }
  
  template<const int lane> 
    __forceinline mic_f lane_shuffle_gather(const mic_f &v0,const mic_f &v1,const mic_f &v2,const mic_f &v3)
  {
    mic_f t = lshuf<lane,lane,lane,lane>(v0);
    t = lshuf<lane,lane,lane,lane>(0xf0,t,v1);
    t = lshuf<lane,lane,lane,lane>(0xf00,t,v2);
    t = lshuf<lane,lane,lane,lane>(0xf000,t,v3);
    return t;
  }
  
  
  __forceinline mic_f mask_min(const mic_m &mask,const mic_f &c, const mic_f &a, const mic_f &b) { 
    return _mm512_mask_gmin_ps(c,mask,a,b); 
  }; 
  
  __forceinline mic_f mask_max(const mic_m &mask,const mic_f &c, const mic_f &a, const mic_f &b) { 
    return _mm512_mask_gmax_ps(c,mask,a,b); 
  }; 
  
  __forceinline mic_f mask_add(const mic_m &mask,mic_f &c, const mic_f &a, const mic_f &b) { 
    return _mm512_mask_add_ps(c,mask,a,b); 
  }; 


  template<int i>
    __forceinline mic_f align_shift_right(const mic_f &a, const mic_f &b)
    {
      return _mm512_castsi512_ps(_mm512_alignr_epi32(_mm512_castps_si512(a),_mm512_castps_si512(b),i)); 
    };

  template<int i>
    __forceinline mic_f mask_align_shift_right(const mic_m &mask,mic_f &c,const mic_f &a, const mic_f &b)
    {
      return _mm512_castsi512_ps(_mm512_mask_alignr_epi32(_mm512_castps_si512(c),mask,_mm512_castps_si512(a),_mm512_castps_si512(b),i)); 
    };

  __forceinline mic_f shl1_zero_extend(const mic_f &a)
  {
    mic_f z = mic_f::zero();
    return mask_align_shift_right<15>(0xfffe,z,a,a);
  }  
  
  __forceinline mic_f prefix_min(const mic_f &a)
  {
    mic_f v = a;
    v = mask_min(0xaaaa,v,v,swizzle(v,_MM_SWIZ_REG_CDAB));
    v = mask_min(0xcccc,v,v,swizzle(v,_MM_SWIZ_REG_BBBB));
    const mic_f shuf_v0 = shuffle(v,(_MM_PERM_ENUM)_MM_SHUF_PERM(2,2,0,0),_MM_SWIZ_REG_DDDD);
    v = mask_min(0xf0f0,v,v,shuf_v0);
    const mic_f shuf_v1 = shuffle(v,(_MM_PERM_ENUM)_MM_SHUF_PERM(1,1,0,0),_MM_SWIZ_REG_DDDD);
    v = mask_min(0xff00,v,v,shuf_v1);
    return v;  
  }
  
  __forceinline mic_f prefix_max(const mic_f &a)
  {
    mic_f v = a;
    v = mask_max(0xaaaa,v,v,swizzle(v,_MM_SWIZ_REG_CDAB));
    v = mask_max(0xcccc,v,v,swizzle(v,_MM_SWIZ_REG_BBBB));
    const mic_f shuf_v0 = shuffle(v,(_MM_PERM_ENUM)_MM_SHUF_PERM(2,2,0,0),_MM_SWIZ_REG_DDDD);
    v = mask_max(0xf0f0,v,v,shuf_v0);
    const mic_f shuf_v1 = shuffle(v,(_MM_PERM_ENUM)_MM_SHUF_PERM(1,1,0,0),_MM_SWIZ_REG_DDDD);
    v = mask_max(0xff00,v,v,shuf_v1);
    return v;  
  }
  
  __forceinline mic_f prefix_sum(const mic_f &a)
  {
    mic_f v = a;
    v = mask_add(0xaaaa,v,v,swizzle(v,_MM_SWIZ_REG_CDAB));
    v = mask_add(0xcccc,v,v,swizzle(v,_MM_SWIZ_REG_BBBB));
    const mic_f shuf_v0 = shuffle(v,(_MM_PERM_ENUM)_MM_SHUF_PERM(2,2,0,0),_MM_SWIZ_REG_DDDD);
    v = mask_add(0xf0f0,v,v,shuf_v0);
    const mic_f shuf_v1 = shuffle(v,(_MM_PERM_ENUM)_MM_SHUF_PERM(1,1,0,0),_MM_SWIZ_REG_DDDD);
    v = mask_add(0xff00,v,v,shuf_v1);
    return v;  
  }
  
  __forceinline mic_f clampz(const mic_f &a, const mic_f &b) 
  { 
    return _min(_max(a,mic_f::zero()),b);
  } 
  
  __forceinline mic_m lz_add(const mic_f &a, const mic_f &b)
  {
    return lz(a+b);
  };

  __forceinline mic_f lz_add(const mic_f &a, const mic_f &b, mic_m &m_res)
  {
    mic_f c = a + b;
    m_res = lz(c);
    return c;
  };
  
  __forceinline mic_m lz_add(const mic_m &m_mask, const mic_f &a, const mic_f &b)
  {
    mic_f v = undefined();
    v = _mm512_mask_add_ps(v,m_mask,a,b);
    return lz(m_mask,v);  
  };
  
  __forceinline mic_f add_sign(const mic_f &a, const mic_f &b, mic_m &m_mask)
  {
    mic_f t = undefined();
    t = _mm512_mask_add_ps(t,m_mask,a,b);
    m_mask = lz(m_mask,t);
    return t;
  }

  __forceinline mic_f _mm512_permutev_ps(__m512i index, mic_f v)
  {
    return _mm512_castsi512_ps(_mm512_permutev_epi32(index,_mm512_castps_si512(v)));  
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Ternary Operators
  ////////////////////////////////////////////////////////////////////////////////

#if defined(__AVX2__)
  __forceinline const mic_f madd  ( const mic_f& a, const mic_f& b, const mic_f& c) { return _mm512_fmadd_ps(a,b,c); }
  __forceinline const mic_f msub  ( const mic_f& a, const mic_f& b, const mic_f& c) { return _mm512_fmsub_ps(a,b,c); }
  __forceinline const mic_f nmadd ( const mic_f& a, const mic_f& b, const mic_f& c) { return _mm512_fnmadd_ps(a,b,c); }
  __forceinline const mic_f nmsub ( const mic_f& a, const mic_f& b, const mic_f& c) { return _mm512_fnmsub_ps(a,b,c); }
#else
  __forceinline const mic_f madd  ( const mic_f& a, const mic_f& b, const mic_f& c) { return a*b+c; }
  __forceinline const mic_f msub  ( const mic_f& a, const mic_f& b, const mic_f& c) { return a*b-c; }
  __forceinline const mic_f nmadd ( const mic_f& a, const mic_f& b, const mic_f& c) { return -a*b-c;}
  __forceinline const mic_f nmsub ( const mic_f& a, const mic_f& b, const mic_f& c) { return c-a*b; }
#endif

  ////////////////////////////////////////////////////////////////////////////////
  /// Select
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const mic_f select( const mic_m& s, const mic_f& t, const mic_f& f ) {
    return _mm512_mask_blend_ps(s, f, t);
  }
  
  __forceinline const mic_f rcp_safe(const mic_f& a) { return rcp(select(a==0.0f,mic_f(1E-10f),a)); }
}

#endif

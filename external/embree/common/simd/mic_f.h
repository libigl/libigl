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

namespace embree
{
  /*! 16-wide MIC float type. */
  class mic_f 
  {
  public:
    
    union  { 
      __m512 v; 
      float f[16]; 
    };

    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////
        
    __forceinline mic_f() {}
    __forceinline mic_f(const mic_f& t) { v = t; };
    __forceinline mic_f& operator=(const mic_f& f) { v = f.v; return *this; };

    __forceinline mic_f(const __m512& t) { v = t; };
    __forceinline operator __m512 () const { return v; };
    
    __forceinline mic_f(const float& f) { 
      v = _mm512_set_1to16_ps(f);
    }
    __forceinline mic_f(const float& a, const float& b, const float& c, const float& d) { 
      v = _mm512_set_4to16_ps(a,b,c,d);  
    }
    
    __forceinline explicit mic_f(const __m512i& a) { 
      v = _mm512_cvtfxpnt_round_adjustepi32_ps(a, _MM_FROUND_NO_EXC,_MM_EXPADJ_NONE);
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline mic_f( ZeroTy   ) : v(_mm512_setzero_ps()) {}
    __forceinline mic_f( OneTy    ) : v(_mm512_set_1to16_ps(1.0f)) {}
    __forceinline mic_f( PosInfTy ) : v(_mm512_set_1to16_ps(pos_inf)) {}
    __forceinline mic_f( NegInfTy ) : v(_mm512_set_1to16_ps(neg_inf)) {}
    __forceinline mic_f( StepTy )   : v(_mm512_set_16to16_ps(15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0)) {}
    __forceinline mic_f( NaNTy    ) : v(_mm512_set_1to16_ps(nan)) {}

    __forceinline static mic_f undefined() { return _mm512_undefined(); }
    __forceinline static mic_f zero() { return _mm512_setzero_ps(); }
    __forceinline static mic_f one () { return _mm512_set_1to16_ps(1.0f); }
    __forceinline static mic_f ulp () { return _mm512_set_1to16_ps(embree::ulp); }
    __forceinline static mic_f inf () { return _mm512_set_1to16_ps((float)pos_inf); }
    __forceinline static mic_f minus_inf () { return _mm512_set_1to16_ps((float)neg_inf); }

    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline float&       operator[](const size_t index)       { return f[index]; };
    __forceinline const float& operator[](const size_t index) const { return f[index]; };
  };

  ////////////////////////////////////////////////////////////////////////////////
  /// Unary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const mic_f cast      (const __m512i& a) { return _mm512_castsi512_ps(a); }
  __forceinline const mic_f operator +( const mic_f& a ) { return a; }
  __forceinline const mic_f operator -( const mic_f& a ) { return _mm512_mul_ps(a,mic_f(-1)); }
  __forceinline const mic_f abs       ( const mic_f& a ) { return _mm512_gmaxabs_ps(a,a); }
  __forceinline const mic_f signmsk   ( const mic_f& a ) { return _mm512_castsi512_ps(_mm512_and_epi32(_mm512_castps_si512(a),_mm512_set1_epi32(0x80000000))); }

  __forceinline const mic_f rcp  ( const mic_f& a ) { return _mm512_rcp23_ps(a); };

  __forceinline const mic_f sqr  ( const mic_f& a ) { return _mm512_mul_ps(a,a); }
  __forceinline const mic_f sqrt ( const mic_f& a ) { return _mm512_sqrt_ps(a); }
  __forceinline const mic_f rsqrt( const mic_f& a ) { return _mm512_invsqrt_ps(a); }
  

  //__forceinline mic_f floor(const mic_f& a) { return _mm512_floor_ps(a); }
  //__forceinline mic_f ceil(const mic_f& a) { return _mm512_ceil_ps(a); }
  __forceinline mic_f trunc(const mic_f& a) { return _mm512_trunc_ps(a); } 
  
  __forceinline mic_f exp(const mic_f& a) { return _mm512_exp_ps(a); }
  __forceinline mic_f exp2(const mic_f& a) { return _mm512_exp2_ps(a); }
  __forceinline mic_f pow(const mic_f& a, mic_f b) { return _mm512_pow_ps(a,b); }
  
  __forceinline mic_f log(const mic_f& a) { return _mm512_log_ps(a); }
  __forceinline mic_f log2(const mic_f& a) { return _mm512_log2_ps(a); }
  __forceinline mic_f log10(const mic_f& a) { return _mm512_log10_ps(a); }
  
  __forceinline mic_f sin(const mic_f& a) { return _mm512_sin_ps(a); } 
  __forceinline mic_f cos(const mic_f& a) { return _mm512_cos_ps(a); }
  __forceinline mic_f tan(const mic_f& a) { return _mm512_tan_ps(a); } 
  
  __forceinline mic_f asin(const mic_f& a) { return _mm512_asin_ps(a); }
  __forceinline mic_f acos(const mic_f& a) { return _mm512_acos_ps(a); }
  __forceinline mic_f atan(const mic_f& a) { return _mm512_atan_ps(a); }
  __forceinline mic_f atan2(const mic_f& a, mic_f b) { return _mm512_atan2_ps(a,b); }
  
  __forceinline mic_f sinh(const mic_f& a) { return _mm512_sinh_ps(a); } 
  __forceinline mic_f cosh(const mic_f& a) { return _mm512_cosh_ps(a); }
  __forceinline mic_f tanh(const mic_f& a) { return _mm512_tan_ps(a); } 
  
  __forceinline mic_f asinh(const mic_f& a) { return _mm512_asinh_ps(a); }
  __forceinline mic_f acosh(const mic_f& a) { return _mm512_acosh_ps(a); }
  __forceinline mic_f atanh(const mic_f& a) { return _mm512_atanh_ps(a); }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const mic_f operator +( const mic_f& a, const mic_f& b ) { return _mm512_add_ps(a, b); }
  __forceinline const mic_f operator +( const mic_f& a, const float& b ) { return a + mic_f(b); }
  __forceinline const mic_f operator +( const float& a, const mic_f& b ) { return mic_f(a) + b; }

  __forceinline const mic_f operator -( const mic_f& a, const mic_f& b ) { return _mm512_sub_ps(a, b); }
  __forceinline const mic_f operator -( const mic_f& a, const float& b ) { return a - mic_f(b); }
  __forceinline const mic_f operator -( const float& a, const mic_f& b ) { return mic_f(a) - b; }

  __forceinline const mic_f operator *( const mic_f& a, const mic_f& b ) { return _mm512_mul_ps(a, b); }
  __forceinline const mic_f operator *( const mic_f& a, const float& b ) { return a * mic_f(b); }
  __forceinline const mic_f operator *( const float& a, const mic_f& b ) { return mic_f(a) * b; }

  __forceinline const mic_f operator /( const mic_f& a, const mic_f& b ) { return _mm512_div_ps(a,b); }
  __forceinline const mic_f operator /( const mic_f& a, const float& b ) { return a/mic_f(b); }
  __forceinline const mic_f operator /( const float& a, const mic_f& b ) { return mic_f(a)/b; }
  
  __forceinline const mic_f operator^(const mic_f& a, const mic_f& b) { 
    return  _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(a),_mm512_castps_si512(b))); 
  }
  
  __forceinline const mic_f min( const mic_f& a, const mic_f& b ) { return _mm512_gmin_ps(a,b); }
  __forceinline const mic_f min( const mic_f& a, const float& b ) { return _mm512_gmin_ps(a,mic_f(b)); }
  __forceinline const mic_f min( const float& a, const mic_f& b ) { return _mm512_gmin_ps(mic_f(a),b); }

  __forceinline const mic_f max( const mic_f& a, const mic_f& b ) { return _mm512_gmax_ps(a,b); }
  __forceinline const mic_f max( const mic_f& a, const float& b ) { return _mm512_gmax_ps(a,mic_f(b)); }
  __forceinline const mic_f max( const float& a, const mic_f& b ) { return _mm512_gmax_ps(mic_f(a),b); }

  __forceinline mic_f mask_add(const mic_m& mask, const mic_f& c, const mic_f& a, const mic_f& b) { return _mm512_mask_add_ps (c,mask,a,b); }; 
  __forceinline mic_f mask_min(const mic_m& mask, const mic_f& c, const mic_f& a, const mic_f& b) { return _mm512_mask_gmin_ps(c,mask,a,b); }; 
  __forceinline mic_f mask_max(const mic_m& mask, const mic_f& c, const mic_f& a, const mic_f& b) { return _mm512_mask_gmax_ps(c,mask,a,b); }; 

  ////////////////////////////////////////////////////////////////////////////////
  /// Ternary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline mic_f madd (const mic_f& a, const mic_f& b, const mic_f& c) { return _mm512_fmadd_ps(a,b,c); }
  __forceinline mic_f msub (const mic_f& a, const mic_f& b, const mic_f& c) { return _mm512_fmsub_ps(a,b,c); }
  __forceinline mic_f nmadd (const mic_f& a, const mic_f& b, const mic_f& c) { return _mm512_fnmadd_ps(a,b,c); }
  __forceinline mic_f nmsub (const mic_f& a, const mic_f& b, const mic_f& c) { return _mm512_fnmsub_ps(a,b,c); }

  __forceinline mic_f mask_msub (const mic_m& mask,const mic_f& a, const mic_f& b, const mic_f& c) { return _mm512_mask_fmsub_ps(a,mask,b,c); }
  
  __forceinline mic_f madd231 (const mic_f& a, const mic_f& b, const mic_f& c) { return _mm512_fmadd_ps(c,b,a); }
  __forceinline mic_f msub213 (const mic_f& a, const mic_f& b, const mic_f& c) { return _mm512_fmsub_ps(a,b,c); }
  __forceinline mic_f msub231 (const mic_f& a, const mic_f& b, const mic_f& c) { return _mm512_fmsub_ps(c,b,a); }
  __forceinline mic_f msubr231(const mic_f& a, const mic_f& b, const mic_f& c) { return _mm512_fnmadd_ps(c,b,a); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline mic_f& operator +=( mic_f& a, const mic_f& b ) { return a = a + b; }
  __forceinline mic_f& operator +=( mic_f& a, const float& b ) { return a = a + b; }
  
  __forceinline mic_f& operator -=( mic_f& a, const mic_f& b ) { return a = a - b; }
  __forceinline mic_f& operator -=( mic_f& a, const float& b ) { return a = a - b; }
  
  __forceinline mic_f& operator *=( mic_f& a, const mic_f& b ) { return a = a * b; }
  __forceinline mic_f& operator *=( mic_f& a, const float& b ) { return a = a * b; }

  __forceinline mic_f& operator /=( mic_f& a, const mic_f& b ) { return a = a / b; }
  __forceinline mic_f& operator /=( mic_f& a, const float& b ) { return a = a / b; }

  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators + Select
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const mic_m operator ==( const mic_f& a, const mic_f& b ) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_EQ); }
  __forceinline const mic_m operator ==( const mic_f& a, const float& b ) { return a == mic_f(b); }
  __forceinline const mic_m operator ==( const float& a, const mic_f& b ) { return mic_f(a) == b; }

  __forceinline const mic_m operator !=( const mic_f& a, const mic_f& b ) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_NE); }
  __forceinline const mic_m operator !=( const mic_f& a, const float& b ) { return a != mic_f(b); }
  __forceinline const mic_m operator !=( const float& a, const mic_f& b ) { return mic_f(a) != b; }

  __forceinline const mic_m operator < ( const mic_f& a, const mic_f& b ) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_LT); }
  __forceinline const mic_m operator < ( const mic_f& a, const float& b ) { return a <  mic_f(b); }
  __forceinline const mic_m operator < ( const float& a, const mic_f& b ) { return mic_f(a) <  b; }

  __forceinline const mic_m operator >=( const mic_f& a, const mic_f& b ) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_GE); }
  __forceinline const mic_m operator >=( const mic_f& a, const float& b ) { return a >= mic_f(b); }
  __forceinline const mic_m operator >=( const float& a, const mic_f& b ) { return mic_f(a) >= b; }

  __forceinline const mic_m operator > ( const mic_f& a, const mic_f& b ) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_GT); }
  __forceinline const mic_m operator > ( const mic_f& a, const float& b ) { return a >  mic_f(b); }
  __forceinline const mic_m operator > ( const float& a, const mic_f& b ) { return mic_f(a) >  b; }

  __forceinline const mic_m operator <=( const mic_f& a, const mic_f& b ) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_LE); }
  __forceinline const mic_m operator <=( const mic_f& a, const float& b ) { return a <= mic_f(b); }
  __forceinline const mic_m operator <=( const float& a, const mic_f& b ) { return mic_f(a) <= b; }

  __forceinline mic_m eq(                   const mic_f& a, const mic_f& b) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_EQ); }
  __forceinline mic_m eq(const mic_m& mask, const mic_f& a, const mic_f& b) { return _mm512_mask_cmp_ps_mask(mask,a,b,_MM_CMPINT_EQ); }

  __forceinline mic_m ne(                   const mic_f& a, const mic_f& b) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_NE); }
  __forceinline mic_m ne(const mic_m& mask, const mic_f& a, const mic_f& b) { return _mm512_mask_cmp_ps_mask(mask,a,b,_MM_CMPINT_NE); }

  __forceinline mic_m lt(                   const mic_f& a, const mic_f& b) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_LT); }
  __forceinline mic_m lt(const mic_m& mask, const mic_f& a, const mic_f& b) { return _mm512_mask_cmp_ps_mask(mask,a,b,_MM_CMPINT_LT); }
 
  __forceinline mic_m ge(                   const mic_f& a, const mic_f& b) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_GE); }
  __forceinline mic_m ge(const mic_m& mask, const mic_f& a, const mic_f& b) { return _mm512_mask_cmp_ps_mask(mask,a,b,_MM_CMPINT_GE); }

  __forceinline mic_m gt(                   const mic_f& a, const mic_f& b) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_GT); }
  __forceinline mic_m gt(const mic_m& mask, const mic_f& a, const mic_f& b) { return _mm512_mask_cmp_ps_mask(mask,a,b,_MM_CMPINT_GT); }
  
  __forceinline mic_m le(                   const mic_f& a, const mic_f& b) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_LE); }
  __forceinline mic_m le(const mic_m& mask, const mic_f& a, const mic_f& b) { return _mm512_mask_cmp_ps_mask(mask,a,b,_MM_CMPINT_LE); }
  
  __forceinline const mic_f select( const mic_m& s, const mic_f& t, const mic_f& f ) {
    return _mm512_mask_blend_ps(s, f, t);
  }


  __forceinline void xchg(mic_m m, mic_f& a, mic_f& b) 
  {
    mic_f c = a;
    a = select(m,b,a);
    b = select(m,c,b); 
  }
  ////////////////////////////////////////////////////////////////////////////////
  /// Rounding Functions
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline mic_f vround(const mic_f& f, 
                             const _MM_ROUND_MODE_ENUM mode, 
                             const _MM_EXP_ADJ_ENUM exp = _MM_EXPADJ_NONE) 
  { 
    return _mm512_round_ps(f,mode,exp); 
  }
  
  __forceinline mic_f floor(const mic_f& a) { return _mm512_round_ps(a,_MM_ROUND_MODE_DOWN, _MM_EXPADJ_NONE); }
  __forceinline mic_f ceil (const mic_f& a) { return _mm512_round_ps(a,_MM_ROUND_MODE_UP  , _MM_EXPADJ_NONE); }

  __forceinline const mic_f rcp_nr  ( const mic_f& a ) { 
    const mic_f ra = _mm512_rcp23_ps(a); 
    return (ra+ra) - (ra * a * ra);
  };

  ////////////////////////////////////////////////////////////////////////////////
  /// Movement/Shifting/Shuffling Functions
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline mic_f swizzle(const mic_f& x,_MM_SWIZZLE_ENUM perm32 ) { return _mm512_swizzle_ps(x,perm32); }
  __forceinline mic_f permute(const mic_f& x,_MM_PERM_ENUM    perm128) { return _mm512_permute4f128_ps(x,perm128); }
  
  template<int D, int C, int B, int A> __forceinline mic_f swizzle   (const mic_f& v) { return _mm512_shuffle_ps(v,_MM_SHUF_PERM(D,C,B,A)); }
  template<int A>                      __forceinline mic_f swizzle   (const mic_f& x) { return swizzle<A,A,A,A>(v); }
  template<>                           __forceinline mic_f swizzle<0>(const mic_f& x) { return swizzle(x,_MM_SWIZ_REG_AAAA); }
  template<>                           __forceinline mic_f swizzle<1>(const mic_f& x) { return swizzle(x,_MM_SWIZ_REG_BBBB); }
  template<>                           __forceinline mic_f swizzle<2>(const mic_f& x) { return swizzle(x,_MM_SWIZ_REG_CCCC); }
  template<>                           __forceinline mic_f swizzle<3>(const mic_f& x) { return swizzle(x,_MM_SWIZ_REG_DDDD); }

  template<int D, int C, int B, int A> __forceinline mic_f permute(const mic_f& v) { return permute(v,_MM_SHUF_PERM(D,C,B,A)); }
  template<int A>                      __forceinline mic_f permute(const mic_f& x) { return permute<A,A,A,A>(x); }

  __forceinline mic_f shuffle(const mic_f& x,_MM_PERM_ENUM perm128, _MM_SWIZZLE_ENUM perm32) { return swizzle(permute(x,perm128),perm32); }
  
  __forceinline mic_f shuffle(const mic_m& mask, mic_f& v, const mic_f& x,_MM_PERM_ENUM perm128, _MM_SWIZZLE_ENUM perm32)  {
    return _mm512_mask_swizzle_ps(_mm512_mask_permute4f128_ps(v,mask,x,perm128),mask,x,perm32);  
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

  __forceinline mic_f _mm512_permutev_ps(__m512i index, mic_f v)
  {
    return _mm512_castsi512_ps(_mm512_permutev_epi32(index,_mm512_castps_si512(v)));  
  }

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


  ////////////////////////////////////////////////////////////////////////////////
  /// Reductions
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline float reduce_add(mic_f a) { return _mm512_reduce_add_ps(a); }
  __forceinline float reduce_mul(mic_f a) { return _mm512_reduce_mul_ps(a); }
  __forceinline float reduce_min(mic_f a) { return _mm512_reduce_min_ps(a); }
  __forceinline float reduce_max(mic_f a) { return _mm512_reduce_max_ps(a); }

  __forceinline mic_f vreduce_min2(mic_f x) {                      return min(x,swizzle(x,_MM_SWIZ_REG_BADC)); }
  __forceinline mic_f vreduce_min4(mic_f x) { x = vreduce_min2(x); return min(x,swizzle(x,_MM_SWIZ_REG_CDAB)); }
  __forceinline mic_f vreduce_min8(mic_f x) { x = vreduce_min4(x); return min(x,permute(x,_MM_SHUF_PERM(2,3,0,1))); }
  __forceinline mic_f vreduce_min (mic_f x) { x = vreduce_min8(x); return min(x,permute(x,_MM_SHUF_PERM(1,0,3,2))); }

  __forceinline mic_f vreduce_max2(mic_f x) {                      return max(x,swizzle(x,_MM_SWIZ_REG_BADC)); }
  __forceinline mic_f vreduce_max4(mic_f x) { x = vreduce_max2(x); return max(x,swizzle(x,_MM_SWIZ_REG_CDAB)); }
  __forceinline mic_f vreduce_max8(mic_f x) { x = vreduce_max4(x); return max(x,permute(x,_MM_SHUF_PERM(2,3,0,1))); }
  __forceinline mic_f vreduce_max (mic_f x) { x = vreduce_max8(x); return max(x,permute(x,_MM_SHUF_PERM(1,0,3,2))); }

  __forceinline mic_f vreduce_add2(mic_f x) {                      return x + swizzle(x,_MM_SWIZ_REG_BADC); }
  __forceinline mic_f vreduce_add4(mic_f x) { x = vreduce_add2(x); return x + swizzle(x,_MM_SWIZ_REG_CDAB); }
  __forceinline mic_f vreduce_add8(mic_f x) { x = vreduce_add4(x); return x + permute(x,_MM_SHUF_PERM(2,3,0,1)); }
  __forceinline mic_f vreduce_add (mic_f x) { x = vreduce_add8(x); return x + permute(x,_MM_SHUF_PERM(1,0,3,2)); }

  __forceinline size_t select_min(const mic_f& v) { return __bsf(movemask(v == vreduce_min(v))); }
  __forceinline size_t select_max(const mic_f& v) { return __bsf(movemask(v == vreduce_max(v))); }

  __forceinline size_t select_min(const mic_m& valid, const mic_f& v, const mic_f &max_value) { const mic_f a = select(valid,v,max_value); return __bsf(movemask(a == vreduce_min(a))); }

  __forceinline size_t select_max(const mic_m& valid, const mic_f& v, const mic_f &min_value) { const mic_f a = select(valid,v,min_value); return __bsf(movemask(a == vreduce_max(a))); }

  __forceinline size_t select_max(const mic_m& valid, const mic_f& v) { const mic_f a = select(valid,v,mic_f(neg_inf)); return __bsf(movemask(valid & (a == vreduce_max(a)))); }

  __forceinline size_t select_min(const mic_m& valid, const mic_f& v) { const mic_f a = select(valid,v,mic_f(pos_inf)); return __bsf(movemask(valid & (a == vreduce_min(a)))); }
  
  __forceinline mic_f prefix_sum(const mic_f& a)
  {
    mic_f v = a;
    v = mask_add(0xaaaa,v,v,swizzle(v,_MM_SWIZ_REG_CDAB));
    v = mask_add(0xcccc,v,v,swizzle(v,_MM_SWIZ_REG_BBBB));
    const mic_f shuf_v0 = shuffle(v,_MM_SHUF_PERM(2,2,0,0),_MM_SWIZ_REG_DDDD);
    v = mask_add(0xf0f0,v,v,shuf_v0);
    const mic_f shuf_v1 = shuffle(v,_MM_SHUF_PERM(1,1,0,0),_MM_SWIZ_REG_DDDD);
    v = mask_add(0xff00,v,v,shuf_v1);
    return v;  
  }

  __forceinline mic_f prefix_min(const mic_f& a)
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
  
  __forceinline mic_f prefix_max(const mic_f& a)
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

  __forceinline mic_f set_min4(mic_f x) {
    x = min(x,swizzle(x,_MM_SWIZ_REG_BADC));
    x = min(x,swizzle(x,_MM_SWIZ_REG_CDAB));
    return x;
  }

  __forceinline mic_f set_min_lanes(mic_f x) {
    x = min(x,_mm512_permute4f128_ps(x, _MM_PERM_CDAB));
    x = min(x,_mm512_permute4f128_ps(x, _MM_PERM_BADC));
    return x;
  }

  __forceinline mic_f set_max_lanes(mic_f x) {
    x = max(x,_mm512_permute4f128_ps(x, _MM_PERM_CDAB));
    x = max(x,_mm512_permute4f128_ps(x, _MM_PERM_BADC));
    return x;
  }

  __forceinline mic_f set_min16(mic_f x) {
    return set_min_lanes(set_min4(x));
  }



  ////////////////////////////////////////////////////////////////////////////////
  /// Memory load and store operations
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline mic_f load1f(const void *f) { 
    return _mm512_extload_ps(f,_MM_UPCONV_PS_NONE,_MM_BROADCAST_1X16,_MM_HINT_NONE);  
  }
  
  __forceinline mic_f load16f(const void *f) { 
    return _mm512_extload_ps(f,_MM_UPCONV_PS_NONE,_MM_BROADCAST_16X16,_MM_HINT_NONE);  
  }

  __forceinline mic_f load16f(const mic_f d, const mic_m valid, const void *f) { 
    return _mm512_mask_extload_ps(d,valid,f,_MM_UPCONV_PS_NONE,_MM_BROADCAST_16X16,_MM_HINT_NONE);  
  }

  __forceinline mic_f broadcast1to16f(const void *f) { 
    return _mm512_extload_ps(f,_MM_UPCONV_PS_NONE,_MM_BROADCAST_1X16,_MM_HINT_NONE);  
  }
  
  __forceinline mic_f broadcast4to16f(const void *f) { 
    return _mm512_extload_ps(f,_MM_UPCONV_PS_NONE,_MM_BROADCAST_4X16,0);  
  }

  __forceinline mic_f broadcast8to16f(const void *f) { 
    return  _mm512_castpd_ps(_mm512_extload_pd(f,_MM_UPCONV_PD_NONE,_MM_BROADCAST_4X8,0));  
  }
    
  __forceinline mic_f load16f_uint8(const unsigned char *const ptr) {
    return _mm512_mul_ps(_mm512_extload_ps(ptr,_MM_UPCONV_PS_UINT8,_MM_BROADCAST_16X16,_MM_HINT_NONE),mic_f(1.0f/255.0f));  
  }

  __forceinline mic_f load16f_uint16(const unsigned short *const ptr) {
    return _mm512_mul_ps(_mm512_extload_ps(ptr,_MM_UPCONV_PS_UINT16,_MM_BROADCAST_16X16,_MM_HINT_NONE),mic_f(1.0f/65535.0f));  
  }

  __forceinline mic_f uload16f_low_uint8(const mic_m& mask, const void* addr, const mic_f& v1) {
    return _mm512_mask_extloadunpacklo_ps(v1, mask, addr, _MM_UPCONV_PS_UINT8, _MM_HINT_NONE);
  }

  __forceinline mic_f load16f_int8(const char *const ptr) {
    return _mm512_mul_ps(_mm512_extload_ps(ptr,_MM_UPCONV_PS_SINT8,_MM_BROADCAST_16X16,_MM_HINT_NONE),mic_f(1.0f/127.0f));  
  }


  __forceinline mic_f gather16f_4f(const float *__restrict__ const ptr0,
                                   const float *__restrict__ const ptr1,
                                   const float *__restrict__ const ptr2,
                                   const float *__restrict__ const ptr3) 
  {
    mic_f v = broadcast4to16f(ptr0);
    v = select((mic_m)0x00f0,broadcast4to16f(ptr1),v);
    v = select((mic_m)0x0f00,broadcast4to16f(ptr2),v);
    v = select((mic_m)0xf000,broadcast4to16f(ptr3),v);
    return v;
  }
  
  __forceinline mic_f gather16f_4f(const mic_f& v0,
                                   const mic_f& v1,
                                   const mic_f& v2,
                                   const mic_f& v3) 
  {
    mic_f v = v0;
    v = select((mic_m)0xf0  ,v1,v);
    v = select((mic_m)0xf00 ,v2,v);
    v = select((mic_m)0xf000,v3,v);
    return v;
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
    
    mic_i v = v_mask &  broadcast4to16i((const int*)ptr0);
    v = mask_and(m_00f0,v,v_mask, broadcast4to16i((const int*)ptr1));
    v = mask_and(m_0f00,v,v_mask, broadcast4to16i((const int*)ptr2));
    v = mask_and(m_f000,v,v_mask, broadcast4to16i((const int*)ptr3));
    return cast(v);
  }

  __forceinline mic_f gather_2f_zlc(const mic_i &v_mask,
				    const mic_m &mask,
                                    const void *__restrict__ const ptr0,
                                    const void *__restrict__ const ptr1) 
  {
    mic_i v = v_mask &  broadcast4to16i((const int*)ptr0);
    v = mask_and(mask,v,v_mask, broadcast4to16i((const int*)ptr1));
    return cast(v);
  }


  __forceinline mic_f gather16f_4f_align(const void *__restrict__ const ptr0,
					 const void *__restrict__ const ptr1,
					 const void *__restrict__ const ptr2,
					 const void *__restrict__ const ptr3) 
  {
    mic_f v = broadcast4to16f(ptr3);
    v = align_shift_right<12>(v,broadcast4to16f(ptr2));
    v = align_shift_right<12>(v,broadcast4to16f(ptr1));
    v = align_shift_right<12>(v,broadcast4to16f(ptr0));
    return v;
  }


  __forceinline mic_f gather16f_4f_align(const mic_f& v0,
					 const mic_f& v1,
					 const mic_f& v2,
					 const mic_f& v3) 
  {
    mic_f v = v3;
    v = align_shift_right<12>(v,v2);
    v = align_shift_right<12>(v,v1);
    v = align_shift_right<12>(v,v0);
    return v;
  }

  __forceinline mic_f gather_4f_zlc_align(const mic_i &v_mask,
					  const void *__restrict__ const ptr0,
					  const void *__restrict__ const ptr1,
					  const void *__restrict__ const ptr2,
					  const void *__restrict__ const ptr3) 
  {
    mic_f v = gather16f_4f_align(ptr0,ptr1,ptr2,ptr3);
    return cast(cast(v) & v_mask);
  }


  __forceinline mic_f uload16f(const mic_m& mask,const float *const addr) {
    mic_f r = mic_f::undefined();
    r =_mm512_mask_extloadunpacklo_ps(r, mask,addr, _MM_UPCONV_PS_NONE, _MM_HINT_NONE);
    r = _mm512_mask_extloadunpackhi_ps(r, mask, addr+16, _MM_UPCONV_PS_NONE, _MM_HINT_NONE);  
    return r;
  }

  __forceinline mic_f uload16f(const float *const addr) {
    mic_f r = mic_f::undefined();
    r =_mm512_extloadunpacklo_ps(r, addr, _MM_UPCONV_PS_NONE, _MM_HINT_NONE);
    return _mm512_extloadunpackhi_ps(r, addr+16, _MM_UPCONV_PS_NONE, _MM_HINT_NONE);  
  }
  
  __forceinline mic_f uload16f_low(const float *const addr) {
    mic_f r = mic_f::undefined();
    return _mm512_extloadunpacklo_ps(r, addr, _MM_UPCONV_PS_NONE, _MM_HINT_NONE);
  }
  
  __forceinline mic_f uload16f_low(const mic_m& mask, const void* addr, const mic_f& v1) {
    return _mm512_mask_extloadunpacklo_ps(v1, mask, addr, _MM_UPCONV_PS_NONE, _MM_HINT_NONE);
  }

  __forceinline mic_f uload16f_low(const mic_m& mask, const void* addr) {
    mic_f v1 = mic_f::undefined();
    return _mm512_mask_extloadunpacklo_ps(v1, mask, addr, _MM_UPCONV_PS_NONE, _MM_HINT_NONE);
  }
  
  __forceinline void ustore16f(float *addr, const mic_f& reg) {
    _mm512_extpackstorelo_ps(addr+0 ,reg, _MM_DOWNCONV_PS_NONE , 0);
    _mm512_extpackstorehi_ps(addr+16 ,reg, _MM_DOWNCONV_PS_NONE , 0);
  }
  
  /* pass by value to avoid compiler generating inefficient code */
  __forceinline void compactustore16f(const mic_m& mask, float *addr, const mic_f reg) {
    _mm512_mask_extpackstorelo_ps(addr+0 ,mask, reg, _MM_DOWNCONV_PS_NONE , 0);
    _mm512_mask_extpackstorehi_ps(addr+16 ,mask, reg, _MM_DOWNCONV_PS_NONE , 0);
  }
  
  __forceinline void compactustore16f_low(const mic_m& mask, float * addr, const mic_f &reg) {
    _mm512_mask_extpackstorelo_ps(addr+0 ,mask, reg, _MM_DOWNCONV_PS_NONE , 0);
  }

  __forceinline void compactustore16f_low_uint8(const mic_m& mask, void * addr, const mic_f &reg) {
    _mm512_mask_extpackstorelo_ps(addr+0 ,mask, reg, _MM_DOWNCONV_PS_UINT8 , 0);
  }
  
  __forceinline void ustore16f_low(float * addr, const mic_f& reg) {
    _mm512_extpackstorelo_ps(addr+0 ,reg, _MM_DOWNCONV_PS_NONE , 0);
  }
  
  __forceinline void compactustore16f_high(const mic_m& mask, float *addr, const mic_f& reg) {
    _mm512_extpackstorehi_ps(addr+0 ,reg, _MM_DOWNCONV_PS_NONE , 0);
  }
  
  __forceinline void store16f(const mic_m& mask, void* addr, const mic_f& v2) {
    _mm512_mask_extstore_ps(addr,mask,v2,_MM_DOWNCONV_PS_NONE,0);
  }
  
  __forceinline void store16f(void* addr, const mic_f& v2) {
    _mm512_extstore_ps(addr,v2,_MM_DOWNCONV_PS_NONE,0);
  }

  __forceinline void store16f_int8(void* addr, const mic_f& v2) {
    _mm512_extstore_ps(addr,v2,_MM_DOWNCONV_PS_SINT8,0);
  }

  __forceinline void store16f_uint16(void* addr, const mic_f& v2) {
    _mm512_extstore_ps(addr,v2,_MM_DOWNCONV_PS_UINT16,0);
  }

  __forceinline void store4f_int8(void* addr, const mic_f& v1) {
    assert((unsigned long)addr % 4 == 0);
    _mm512_mask_extpackstorelo_ps(addr,0xf, v1, _MM_DOWNCONV_PS_SINT8 , 0);
  }
  
  __forceinline void store4f(void* addr, const mic_f& v1) {
    assert((unsigned long)addr % 16 == 0);
    _mm512_mask_extpackstorelo_ps(addr,0xf, v1, _MM_DOWNCONV_PS_NONE , 0);
  }

  __forceinline void store3f(void* addr, const mic_f& v1) {
    assert((unsigned long)addr % 16 == 0);
    _mm512_mask_extpackstorelo_ps(addr,0x7, v1, _MM_DOWNCONV_PS_NONE , 0);
  }

  __forceinline void store4f_nt(void* addr, const mic_f& v1) {
    assert((unsigned long)addr % 16 == 0);
    _mm512_mask_extpackstorelo_ps(addr,0xf, v1, _MM_DOWNCONV_PS_NONE , _MM_HINT_NT);
  }
  
  __forceinline void store16f_nt(void *__restrict__ ptr, const mic_f& a) {
    _mm512_storenr_ps(ptr,a);
  }
  
  __forceinline void store16f_ngo(void *__restrict__ ptr, const mic_f& a) {
    _mm512_storenrngo_ps(ptr,a);
  }

  __forceinline void store1f(void *addr, const mic_f& reg) {
    _mm512_mask_extpackstorelo_ps((float*)addr+0  ,(mic_m)1, reg, _MM_DOWNCONV_PS_NONE, _MM_HINT_NONE);
  }
  
  __forceinline mic_f gather16f(const mic_m& mask, const float *const ptr, __m512i index, const _MM_INDEX_SCALE_ENUM scale) {
    mic_f r = mic_f::undefined();
    return _mm512_mask_i32extgather_ps(r,mask,index,ptr,_MM_UPCONV_PS_NONE,scale,0);
  }
  
  __forceinline void scatter16f(const mic_m& mask,const float *const ptr, const __m512i index,const mic_f v, const _MM_INDEX_SCALE_ENUM scale) {
    _mm512_mask_i32extscatter_ps((void*)ptr,mask,index,v,_MM_DOWNCONV_PS_NONE,scale,0);
  }

  __forceinline mic_f loadAOS4to16f(const float& x,const float& y, const float& z)
  {
    mic_f f = mic_f::zero();
    f = select(0x1111,broadcast1to16f(&x),f);
    f = select(0x2222,broadcast1to16f(&y),f);
    f = select(0x4444,broadcast1to16f(&z),f);
    return f;
  }

  __forceinline mic_f loadAOS4to16f(const unsigned int index,
				    const mic_f &x,
				    const mic_f &y,
				    const mic_f &z)
  {
    mic_f f = mic_f::zero();
    f = select(0x1111,broadcast1to16f((float*)&x + index),f);
    f = select(0x2222,broadcast1to16f((float*)&y + index),f);
    f = select(0x4444,broadcast1to16f((float*)&z + index),f);
    return f;
  }

  __forceinline mic_f loadAOS4to16f(const unsigned int index,
				    const mic_f &x,
				    const mic_f &y,
				    const mic_f &z,
				    const mic_f &fill)
  {
    mic_f f = fill;
    f = select(0x1111,broadcast1to16f((float*)&x + index),f);
    f = select(0x2222,broadcast1to16f((float*)&y + index),f);
    f = select(0x4444,broadcast1to16f((float*)&z + index),f);
    return f;
  }

  __forceinline mic_f rcp_safe( const mic_f& a ) { return select(a != mic_f::zero(),_mm512_rcp23_ps(a),mic_f(1E-10f)); };

  ////////////////////////////////////////////////////////////////////////////////
  /// Euclidian Space Operators
  ////////////////////////////////////////////////////////////////////////////////

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
  
  __forceinline mic_f ldot3_zxy(const mic_f &ao, const mic_f &normal) 
  {
    mic_f vv = ao * swizzle(normal,_MM_SWIZ_REG_DACB);
    vv = _mm512_add_ps(vv,swizzle(vv,_MM_SWIZ_REG_CDAB));
    vv = _mm512_add_ps(vv,swizzle(vv,_MM_SWIZ_REG_BADC));
    return vv;        
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Output Operators
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline std::ostream &operator<<(std::ostream& cout, const mic_f& v)
  {
    cout << "<" << v[0];
    for (int i=1; i<16; i++) cout << ", " << v[i];
    cout << ">";
    return cout;
  }
}

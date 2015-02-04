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
  /*! 16-wide MIC integer type. */
  class mic_i
  {
  public:
    
    union  { 
      __m512i v; 
      int32 i[16]; 
    };
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////
       
    __forceinline mic_i() {};
    __forceinline mic_i(const mic_i& t) { v = t.v; };
    __forceinline mic_i& operator=(const mic_i& f) { v = f.v; return *this; };

    __forceinline mic_i(const __m512i& t) { v = t; };
    __forceinline operator __m512i () const { return v; };

    __forceinline mic_i(const int& i) { 
      v = _mm512_set_1to16_epi32(i);
    }
    
    __forceinline mic_i(const int& a, const int& b, const int& c, const int& d) { 
      v = _mm512_set_4to16_epi32(a,b,c,d);      
    }
   
    __forceinline explicit mic_i(const __m512 f) { 
      v = _mm512_cvtfxpnt_round_adjustps_epi32(f,_MM_FROUND_FLOOR,_MM_EXPADJ_NONE);
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline mic_i( ZeroTy   ) : v(_mm512_setzero_epi32()) {}
    __forceinline mic_i( OneTy    ) : v(_mm512_set_1to16_epi32(1)) {}
    __forceinline mic_i( PosInfTy ) : v(_mm512_set_1to16_epi32(pos_inf)) {}
    __forceinline mic_i( NegInfTy ) : v(_mm512_set_1to16_epi32(neg_inf)) {}
    __forceinline mic_i( StepTy )   : v(_mm512_set_16to16_epi32(15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0)) {}

    __forceinline static mic_i zero() { return _mm512_setzero_epi32(); }
    __forceinline static mic_i one () { return _mm512_set_1to16_epi32(1); }
    __forceinline static mic_i neg_one () { return _mm512_set_1to16_epi32(-1); }

    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline int&       operator[](const size_t index)       { return i[index]; };
    __forceinline const int& operator[](const size_t index) const { return i[index]; };

    __forceinline unsigned int&       uint(const size_t index) const      { assert(index < 16); return ((unsigned int*)i)[index]; };
    __forceinline size_t&             uint64(const size_t index)  const     { assert(index < 8); return ((size_t*)i)[index]; };


  };
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Unary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const mic_i cast      ( const __m512& a) { return _mm512_castps_si512(a); }
  __forceinline const mic_i operator +( const mic_i& a ) { return a; }
  __forceinline const mic_i operator -( const mic_i& a ) { return _mm512_sub_epi32(_mm512_setzero_epi32(), a); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const mic_i operator +( const mic_i& a, const mic_i& b ) { return _mm512_add_epi32(a, b); }
  __forceinline const mic_i operator +( const mic_i& a, const int32&  b ) { return a + mic_i(b); }
  __forceinline const mic_i operator +( const int32& a, const mic_i& b ) { return mic_i(a) + b; }

  __forceinline const mic_i operator -( const mic_i& a, const mic_i& b ) { return _mm512_sub_epi32(a, b); }
  __forceinline const mic_i operator -( const mic_i& a, const int32&  b ) { return a - mic_i(b); }
  __forceinline const mic_i operator -( const int32& a, const mic_i& b ) { return mic_i(a) - b; }

  __forceinline const mic_i operator *( const mic_i& a, const mic_i& b ) { return _mm512_mullo_epi32(a, b); }
  __forceinline const mic_i operator *( const mic_i& a, const int32&  b ) { return a * mic_i(b); }
  __forceinline const mic_i operator *( const int32& a, const mic_i& b ) { return mic_i(a) * b; }

  __forceinline const mic_i operator &( const mic_i& a, const mic_i& b ) { return _mm512_and_epi32(a, b); }
  __forceinline const mic_i operator &( const mic_i& a, const int32&  b ) { return a & mic_i(b); }
  __forceinline const mic_i operator &( const int32& a, const mic_i& b ) { return mic_i(a) & b; }

  __forceinline const mic_i operator |( const mic_i& a, const mic_i& b ) { return _mm512_or_epi32(a, b); }
  __forceinline const mic_i operator |( const mic_i& a, const int32&  b ) { return a | mic_i(b); }
  __forceinline const mic_i operator |( const int32& a, const mic_i& b ) { return mic_i(a) | b; }

  __forceinline const mic_i operator ^( const mic_i& a, const mic_i& b ) { return _mm512_xor_epi32(a, b); }
  __forceinline const mic_i operator ^( const mic_i& a, const int32&  b ) { return a ^ mic_i(b); }
  __forceinline const mic_i operator ^( const int32& a, const mic_i& b ) { return mic_i(a) ^ b; }

  __forceinline const mic_i operator <<( const mic_i& a, const int32& n ) { return _mm512_slli_epi32(a, n); }
  __forceinline const mic_i operator >>( const mic_i& a, const int32& n ) { return _mm512_srai_epi32(a, n); }

  __forceinline const mic_i operator <<( const mic_i& a, const mic_i& n ) { return _mm512_sllv_epi32(a, n); }
  __forceinline const mic_i operator >>( const mic_i& a, const mic_i& n ) { return _mm512_srav_epi32(a, n); }

  __forceinline const mic_i sra ( const mic_i& a, const int32& b ) { return _mm512_srai_epi32(a, b); }
  __forceinline const mic_i srl ( const mic_i& a, const int32& b ) { return _mm512_srli_epi32(a, b); }
  
  __forceinline const mic_i min( const mic_i& a, const mic_i& b ) { return _mm512_min_epi32(a, b); }
  __forceinline const mic_i min( const mic_i& a, const int32&  b ) { return min(a,mic_i(b)); }
  __forceinline const mic_i min( const int32&  a, const mic_i& b ) { return min(mic_i(a),b); }

  __forceinline const mic_i max( const mic_i& a, const mic_i& b ) { return _mm512_max_epi32(a, b); }
  __forceinline const mic_i max( const mic_i& a, const int32&  b ) { return max(a,mic_i(b)); }
  __forceinline const mic_i max( const int32&  a, const mic_i& b ) { return max(mic_i(a),b); }
  
  __forceinline const mic_i mask_add(const mic_m& mask, mic_i& c, const mic_i& a, const mic_i& b) { return _mm512_mask_add_epi32(c,mask,a,b); }; 
  __forceinline const mic_i mask_sub(const mic_m& mask, mic_i& c, const mic_i& a, const mic_i& b) { return _mm512_mask_sub_epi32(c,mask,a,b); }; 

  __forceinline const mic_i mask_and(const mic_m& m,mic_i& c, const mic_i& a, const mic_i& b) { return _mm512_mask_and_epi32(c,m,a,b); };
  __forceinline const mic_i mask_or (const mic_m& m,mic_i& c, const mic_i& a, const mic_i& b) { return _mm512_mask_or_epi32(c,m,a,b); };
 
  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline mic_i& operator +=( mic_i& a, const mic_i& b ) { return a = a + b; }
  __forceinline mic_i& operator +=( mic_i& a, const int32&  b ) { return a = a + b; }
  
  __forceinline mic_i& operator -=( mic_i& a, const mic_i& b ) { return a = a - b; }
  __forceinline mic_i& operator -=( mic_i& a, const int32&  b ) { return a = a - b; }

  __forceinline mic_i& operator *=( mic_i& a, const mic_i& b ) { return a = a * b; }
  __forceinline mic_i& operator *=( mic_i& a, const int32&  b ) { return a = a * b; }
  
  __forceinline mic_i& operator &=( mic_i& a, const mic_i& b ) { return a = a & b; }
  __forceinline mic_i& operator &=( mic_i& a, const int32&  b ) { return a = a & b; }
  
  __forceinline mic_i& operator |=( mic_i& a, const mic_i& b ) { return a = a | b; }
  __forceinline mic_i& operator |=( mic_i& a, const int32&  b ) { return a = a | b; }
  
  __forceinline mic_i& operator <<=( mic_i& a, const int32&  b ) { return a = a << b; }
  __forceinline mic_i& operator >>=( mic_i& a, const int32&  b ) { return a = a >> b; }


  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators + Select
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const mic_m operator ==( const mic_i& a, const mic_i& b ) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_EQ); }
  __forceinline const mic_m operator ==( const mic_i& a, const int32& b ) { return a == mic_i(b); }
  __forceinline const mic_m operator ==( const int32& a, const mic_i& b ) { return mic_i(a) == b; }
  
  __forceinline const mic_m operator !=( const mic_i& a, const mic_i& b ) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_NE); }
  __forceinline const mic_m operator !=( const mic_i& a, const int32& b ) { return a != mic_i(b); }
  __forceinline const mic_m operator !=( const int32& a, const mic_i& b ) { return mic_i(a) != b; }
  
  __forceinline const mic_m operator < ( const mic_i& a, const mic_i& b ) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_LT); }
  __forceinline const mic_m operator < ( const mic_i& a, const int32& b ) { return a <  mic_i(b); }
  __forceinline const mic_m operator < ( const int32& a, const mic_i& b ) { return mic_i(a) <  b; }
  
  __forceinline const mic_m operator >=( const mic_i& a, const mic_i& b ) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_GE); }
  __forceinline const mic_m operator >=( const mic_i& a, const int32& b ) { return a >= mic_i(b); }
  __forceinline const mic_m operator >=( const int32& a, const mic_i& b ) { return mic_i(a) >= b; }

  __forceinline const mic_m operator > ( const mic_i& a, const mic_i& b ) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_GT); }
  __forceinline const mic_m operator > ( const mic_i& a, const int32& b ) { return a >  mic_i(b); }
  __forceinline const mic_m operator > ( const int32& a, const mic_i& b ) { return mic_i(a) >  b; }

  __forceinline const mic_m operator <=( const mic_i& a, const mic_i& b ) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_LE); }
  __forceinline const mic_m operator <=( const mic_i& a, const int32& b ) { return a <= mic_i(b); }
  __forceinline const mic_m operator <=( const int32& a, const mic_i& b ) { return mic_i(a) <= b; }

  __forceinline mic_m eq(                  const mic_i& a, const mic_i& b) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_EQ); };
  __forceinline mic_m eq(const mic_m mask, const mic_i& a, const mic_i& b) { return _mm512_mask_cmp_epi32_mask(mask,a,b,_MM_CMPINT_EQ);  };
  
  __forceinline mic_m ne(                  const mic_i& a, const mic_i& b) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_NE); };
  __forceinline mic_m ne(const mic_m mask, const mic_i& a, const mic_i& b) { return _mm512_mask_cmp_epi32_mask(mask,a,b,_MM_CMPINT_NE); };

  __forceinline mic_m lt(                  const mic_i& a, const mic_i& b) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_LT); };
  __forceinline mic_m lt(const mic_m mask, const mic_i& a, const mic_i& b) { return _mm512_mask_cmp_epi32_mask(mask,a,b,_MM_CMPINT_LT); };
 
  __forceinline mic_m ge(                  const mic_i& a, const mic_i& b) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_GE); };
  __forceinline mic_m ge(const mic_m mask, const mic_i& a, const mic_i& b) { return _mm512_mask_cmp_epi32_mask(mask,a,b,_MM_CMPINT_GE); };
  
  __forceinline mic_m gt(                  const mic_i& a, const mic_i& b) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_GT); };
  __forceinline mic_m gt(const mic_m mask, const mic_i& a, const mic_i& b) { return _mm512_mask_cmp_epi32_mask(mask,a,b,_MM_CMPINT_GT); };
  
  __forceinline mic_m le(                  const mic_i& a, const mic_i& b) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_LE); };
  __forceinline mic_m le(const mic_m mask, const mic_i& a, const mic_i& b) { return _mm512_mask_cmp_epi32_mask(mask,a,b,_MM_CMPINT_LE); };
    
 
  __forceinline const mic_i select( const mic_m& m, const mic_i& t, const mic_i& f ) { 
    return _mm512_mask_or_epi32(f,m,t,t); 
  }

  __forceinline void xchg(const mic_m &m, mic_i& a, mic_i& b) { 
    const mic_i c = a; a = select(m,b,a); b = select(m,c,b);  
  }

  __forceinline mic_m test(const mic_m &m, const mic_i& a, const mic_i& b) { 
    return _mm512_mask_test_epi32_mask(m,a,b);
  }

  __forceinline mic_m test(const mic_i& a, const mic_i& b) { 
    return _mm512_test_epi32_mask(a,b);
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Movement/Shifting/Shuffling Functions
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline mic_i swizzle(const mic_i& x,_MM_SWIZZLE_ENUM perm32 ) { return _mm512_swizzle_epi32(x,perm32); }
  __forceinline mic_i permute(const mic_i& x,_MM_PERM_ENUM    perm128) { return _mm512_permute4f128_epi32(x,perm128); }
  
  template<int D, int C, int B, int A> __forceinline mic_i swizzle   (const mic_i& v) { return _mm512_shuffle_epi32(v,_MM_SHUF_PERM(D,C,B,A)); }
  template<int A>                      __forceinline mic_i swizzle   (const mic_i& x) { return swizzle<A,A,A,A>(v); }
  template<>                           __forceinline mic_i swizzle<0>(const mic_i& x) { return swizzle(x,_MM_SWIZ_REG_AAAA); }
  template<>                           __forceinline mic_i swizzle<1>(const mic_i& x) { return swizzle(x,_MM_SWIZ_REG_BBBB); }
  template<>                           __forceinline mic_i swizzle<2>(const mic_i& x) { return swizzle(x,_MM_SWIZ_REG_CCCC); }
  template<>                           __forceinline mic_i swizzle<3>(const mic_i& x) { return swizzle(x,_MM_SWIZ_REG_DDDD); }

  template<int D, int C, int B, int A> __forceinline mic_i permute(const mic_i& v) { return permute(v,_MM_SHUF_PERM(D,C,B,A)); }
  template<int A>                      __forceinline mic_i permute(const mic_i& x) { return permute<A,A,A,A>(x); }

  __forceinline mic_i shuffle(const mic_i& x,_MM_PERM_ENUM    perm128, _MM_SWIZZLE_ENUM perm32) { return swizzle(permute(x,perm128),perm32); }
  
  __forceinline mic_i shuffle(const mic_m& mask, mic_i& v, const mic_i& x,_MM_PERM_ENUM perm128, _MM_SWIZZLE_ENUM perm32)  {
    return _mm512_mask_swizzle_epi32(_mm512_mask_permute4f128_epi32(v,mask,x,perm128),mask,x,perm32);  
  }

  __forceinline mic_i swAAAA(const mic_i &x) {
    return swizzle(x,_MM_SWIZ_REG_AAAA);
  }

  __forceinline mic_i swBBBB(const mic_i &x) {
    return swizzle(x,_MM_SWIZ_REG_BBBB);
  }

  __forceinline mic_i swCCCC(const mic_i &x) {
    return swizzle(x,_MM_SWIZ_REG_CCCC);
  }

  __forceinline mic_i swDDDD(const mic_i &x) {
    return swizzle(x,_MM_SWIZ_REG_DDDD);
  }

  template<int i>
  __forceinline mic_i align_shift_right(const mic_i &a, const mic_i &b)
  {
    return _mm512_alignr_epi32(a,b,i); 
  };

  
  ////////////////////////////////////////////////////////////////////////////////
  /// Reductions
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline int reduce_add(mic_i a) { return _mm512_reduce_add_epi32(a); }
  __forceinline int reduce_mul(mic_i a) { return _mm512_reduce_mul_epi32(a); }
  __forceinline int reduce_min(mic_i a) { return _mm512_reduce_min_epi32(a); }
  __forceinline int reduce_max(mic_i a) { return _mm512_reduce_max_epi32(a); }
  __forceinline int reduce_and(mic_i a) { return _mm512_reduce_and_epi32(a); }
  
  __forceinline mic_i vreduce_min2(mic_i x) {                      return min(x,swizzle(x,_MM_SWIZ_REG_BADC)); }
  __forceinline mic_i vreduce_min4(mic_i x) { x = vreduce_min2(x); return min(x,swizzle(x,_MM_SWIZ_REG_CDAB)); }
  __forceinline mic_i vreduce_min8(mic_i x) { x = vreduce_min4(x); return min(x,permute(x,_MM_SHUF_PERM(2,3,0,1))); }
  __forceinline mic_i vreduce_min (mic_i x) { x = vreduce_min8(x); return min(x,permute(x,_MM_SHUF_PERM(1,0,3,2))); }

  __forceinline mic_i vreduce_max2(mic_i x) {                      return max(x,swizzle(x,_MM_SWIZ_REG_BADC)); }
  __forceinline mic_i vreduce_max4(mic_i x) { x = vreduce_max2(x); return max(x,swizzle(x,_MM_SWIZ_REG_CDAB)); }
  __forceinline mic_i vreduce_max8(mic_i x) { x = vreduce_max4(x); return max(x,permute(x,_MM_SHUF_PERM(2,3,0,1))); }
  __forceinline mic_i vreduce_max (mic_i x) { x = vreduce_max8(x); return max(x,permute(x,_MM_SHUF_PERM(1,0,3,2))); }

  __forceinline mic_i vreduce_and2(mic_i x) {                      return x & swizzle(x,_MM_SWIZ_REG_BADC); }
  __forceinline mic_i vreduce_and4(mic_i x) { x = vreduce_and2(x); return x & swizzle(x,_MM_SWIZ_REG_CDAB); }
  __forceinline mic_i vreduce_and8(mic_i x) { x = vreduce_and4(x); return x & permute(x,_MM_SHUF_PERM(2,3,0,1)); }
  __forceinline mic_i vreduce_and (mic_i x) { x = vreduce_and8(x); return x & permute(x,_MM_SHUF_PERM(1,0,3,2)); }

  __forceinline mic_i vreduce_or2(mic_i x) {                     return x | swizzle(x,_MM_SWIZ_REG_BADC); }
  __forceinline mic_i vreduce_or4(mic_i x) { x = vreduce_or2(x); return x | swizzle(x,_MM_SWIZ_REG_CDAB); }
  __forceinline mic_i vreduce_or8(mic_i x) { x = vreduce_or4(x); return x | permute(x,_MM_SHUF_PERM(2,3,0,1)); }
  __forceinline mic_i vreduce_or (mic_i x) { x = vreduce_or8(x); return x | permute(x,_MM_SHUF_PERM(1,0,3,2)); }

  __forceinline mic_i vreduce_add2(mic_i x) {                      return x + swizzle(x,_MM_SWIZ_REG_BADC); }
  __forceinline mic_i vreduce_add4(mic_i x) { x = vreduce_add2(x); return x + swizzle(x,_MM_SWIZ_REG_CDAB); }
  __forceinline mic_i vreduce_add8(mic_i x) { x = vreduce_add4(x); return x + permute(x,_MM_SHUF_PERM(2,3,0,1)); }
  __forceinline mic_i vreduce_add (mic_i x) { x = vreduce_add8(x); return x + permute(x,_MM_SHUF_PERM(1,0,3,2)); }
  
  __forceinline mic_i prefix_sum(const mic_i& a)
  {
    mic_i v = a;
    v = mask_add(0xaaaa,v,v,swizzle(v,_MM_SWIZ_REG_CDAB));
    v = mask_add(0xcccc,v,v,swizzle(v,_MM_SWIZ_REG_BBBB));
    const mic_i shuf_v0 = shuffle(v,(_MM_PERM_ENUM)_MM_SHUF_PERM(2,2,0,0),_MM_SWIZ_REG_DDDD);
    v = mask_add(0xf0f0,v,v,shuf_v0);
    const mic_i shuf_v1 = shuffle(v,(_MM_PERM_ENUM)_MM_SHUF_PERM(1,1,0,0),_MM_SWIZ_REG_DDDD);
    v = mask_add(0xff00,v,v,shuf_v1);
    return v;  
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Memory load and store operations
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline mic_i load1i(const int *const ptr) { 
    return _mm512_extload_epi32(ptr,_MM_UPCONV_EPI32_NONE,_MM_BROADCAST_1X16,_MM_HINT_NONE);
  }

  __forceinline mic_i load16i(const int *const ptr) {
    return _mm512_extload_epi32(ptr,_MM_UPCONV_EPI32_NONE,_MM_BROADCAST_16X16,_MM_HINT_NONE);
  }
  
  __forceinline mic_i load1i_uint8(const unsigned char *const ptr) { 
    return _mm512_extload_epi32(ptr,_MM_UPCONV_EPI32_UINT8,_MM_BROADCAST_1X16,_MM_HINT_NONE);
  }

  __forceinline mic_i load16i_uint8(const unsigned char *const ptr) {
    return _mm512_extload_epi32(ptr,_MM_UPCONV_EPI32_UINT8,_MM_BROADCAST32_NONE,_MM_HINT_NONE);
  }
  
  __forceinline mic_i broadcast4to16i(const int *const ptr) {
    return _mm512_extload_epi32(ptr,_MM_UPCONV_EPI32_NONE,_MM_BROADCAST_4X16,_MM_HINT_NONE);
  }

  __forceinline mic_i broadcast1to16i(const int *const ptr) {
    return _mm512_extload_epi32(ptr,_MM_UPCONV_EPI32_NONE,_MM_BROADCAST_1X16,_MM_HINT_NONE);
  }  
  
  __forceinline mic_i uload16i(const int *const addr) {
    mic_i r = _mm512_undefined_epi32();
    r =_mm512_extloadunpacklo_epi32(r, addr, _MM_UPCONV_EPI32_NONE, _MM_HINT_NONE);
    return _mm512_extloadunpackhi_epi32(r, addr+16, _MM_UPCONV_EPI32_NONE, _MM_HINT_NONE);  
  }
  
  __forceinline mic_i uload16i_low(const mic_m& mask, const void* addr) {
    mic_i v = _mm512_undefined_epi32();
    return _mm512_mask_extloadunpacklo_epi32(v, mask, addr, _MM_UPCONV_EPI32_NONE, _MM_HINT_NONE);
  }

  __forceinline mic_i uload16i(const mic_m& mask,const int *const addr) {
    mic_i r = _mm512_undefined_epi32();
    r =_mm512_mask_extloadunpacklo_epi32(r, mask,addr, _MM_UPCONV_EPI32_NONE, _MM_HINT_NONE);
    return _mm512_mask_extloadunpackhi_epi32(r, mask,addr+16, _MM_UPCONV_EPI32_NONE, _MM_HINT_NONE);  
  }
  

  __forceinline mic_i gather16i_4i(const int *__restrict__ const ptr0,
                                   const int *__restrict__ const ptr1,
                                   const int *__restrict__ const ptr2,
                                   const int *__restrict__ const ptr3) 
  {
    mic_i v =  broadcast4to16i(ptr0);
    v = select((mic_m)0xf0  , broadcast4to16i(ptr1),v);
    v = select((mic_m)0xf00 , broadcast4to16i(ptr2),v);
    v = select((mic_m)0xf000, broadcast4to16i(ptr3),v);
    return v;
  }


  __forceinline mic_i gather16i_4i_align(const void *__restrict__ const ptr0,
					 const void *__restrict__ const ptr1,
					 const void *__restrict__ const ptr2,
					 const void *__restrict__ const ptr3) 
  {
    mic_i v = broadcast4to16i((const int*)ptr3);
    v = align_shift_right<12>(v,broadcast4to16i((const int*)ptr2));
    v = align_shift_right<12>(v,broadcast4to16i((const int*)ptr1));
    v = align_shift_right<12>(v,broadcast4to16i((const int*)ptr0));
    return v;
  }

  __forceinline mic_i gather16i_4i_align(const mic_i &v0,
					 const mic_i &v1,
					 const mic_i &v2,
					 const mic_i &v3)
  {
    mic_i v = v3;
    v = align_shift_right<12>(v,v2);
    v = align_shift_right<12>(v,v1);
    v = align_shift_right<12>(v,v0);
    return v;
  }

  
  __forceinline mic_i gather16i(const mic_m& mask, const int *const ptr, const mic_i& index,const _MM_INDEX_SCALE_ENUM scale) {
    return _mm512_mask_i32extgather_epi32(_mm512_undefined_epi32(),mask,index,ptr,_MM_UPCONV_EPI32_NONE,scale,0);
  }
  
  __forceinline mic_i gather16i(const mic_m& mask, mic_i& dest, const int *const ptr, const mic_i& index,const _MM_INDEX_SCALE_ENUM scale) {
    return _mm512_mask_i32extgather_epi32(dest,mask,index,ptr,_MM_UPCONV_EPI32_NONE,scale,0);
  }
  
  __forceinline void scatter16i(const mic_m& mask,int *const ptr, const mic_i& index,const mic_i& v, const _MM_INDEX_SCALE_ENUM scale) {
    _mm512_mask_i32extscatter_epi32((int*)ptr,mask,index,v,_MM_DOWNCONV_EPI32_NONE,scale,0);
  }
  
  __forceinline void ustore16i(void *addr, const mic_i& reg) {
    _mm512_extpackstorelo_epi32((int*)addr+0  ,reg, _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
    _mm512_extpackstorehi_epi32((int*)addr+16 ,reg, _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
  }
  
  /* pass by value to avoid compiler generating inefficient code */
  __forceinline void compactustore16i(const mic_m mask,void * addr, const mic_i reg) {
    _mm512_mask_extpackstorelo_epi32((int*)addr+0  ,mask, reg, _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
    _mm512_mask_extpackstorehi_epi32((int*)addr+16 ,mask, reg, _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
  }


  __forceinline void compactustore16i_low(const mic_m mask, void *addr, const mic_i &reg) {
    _mm512_mask_extpackstorelo_epi32((int*)addr+0  ,mask, reg, _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
  }

  __forceinline void store1i(void *addr, const mic_i& reg) {
    _mm512_mask_extpackstorelo_epi32((int*)addr+0  ,(mic_m)1, reg, _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
  }
  
  __forceinline void ustore16i_low(void *addr, const mic_i& reg) {
    _mm512_extpackstorelo_epi32((int*)addr+0  ,reg, _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
  }
  
  __forceinline void compactustore16i_high(const mic_m mask, int *addr, const mic_i& reg) {
    _mm512_mask_extpackstorehi_epi32((int*)addr+16  ,mask, reg, _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
  }
  
  __forceinline void store16i_nr(void *__restrict__ ptr, const mic_i& a) {
    _mm512_storenr_ps(ptr,_mm512_castsi512_ps(a));
  }
  
  __forceinline void store16i_ngo(void *__restrict__ ptr, const mic_i& a) {
    _mm512_storenrngo_ps(ptr,_mm512_castsi512_ps(a));
  }
  
  __forceinline void store16i(const mic_m& mask, void* __restrict__ addr, const mic_i& v2) {
    _mm512_mask_extstore_epi32(addr,mask,v2,_MM_DOWNCONV_EPI32_NONE,_MM_HINT_NONE);
  }

  __forceinline void store16i_nt(void* __restrict__ addr, const mic_i& v2) {
    _mm512_extstore_epi32(addr,v2,_MM_DOWNCONV_EPI32_NONE,_MM_HINT_NT);
  }
  
  __forceinline void store16i(void* __restrict__ addr, const mic_i& v2) {
    _mm512_extstore_epi32(addr,v2,_MM_DOWNCONV_EPI32_NONE,_MM_HINT_NONE);
  }
  
  __forceinline void store16i_uint8(const mic_m& mask, void* __restrict__ addr, const mic_i& v2) {
    _mm512_mask_extstore_epi32(addr,mask,v2,_MM_DOWNCONV_EPI32_UINT8,_MM_HINT_NONE);
  }

  __forceinline mic_i convert_uint32(const __m512 f) { 
    return _mm512_cvtfxpnt_round_adjustps_epu32(f,_MM_FROUND_TO_ZERO,_MM_EXPADJ_NONE);
  }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Output Operators
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline std::ostream &operator<<(std::ostream& cout, const mic_i& v)
    {
      cout << "<" << v[0];
      for (int i=1; i<16; i++) cout << ", " << v[i];
      cout << ">";
      return cout;
    }
}

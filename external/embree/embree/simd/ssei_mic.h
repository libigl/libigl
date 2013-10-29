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

#ifndef __EMBREE_SSEI_MIC_H__
#define __EMBREE_SSEI_MIC_H__

namespace embree
{
  /* memory representation as 4 aligned ints */
  struct __align(16) ssei_m 
  {
    typedef sseb_m Mask;
    typedef ssei_m Int;
    typedef ssef_m Float;
    enum { size = 4 };  // number of SIMD elements
    int i[4]; 
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////
  
    __forceinline ssei_m ( ) {}
    __forceinline ssei_m ( int a                     ) { i[0] = a;  i[1] = a; i[2] = a; i[3] = a; }
    __forceinline ssei_m ( int a, int b, int c, int d) { i[0] = a;  i[1] = b; i[2] = c; i[3] = d; }

    __forceinline ssei_m ( const __m512i& other ) { 
      assert((size_t)this % 16 == 0); 
      _mm512_mask_extpackstorelo_epi32(this,0xf,other,_MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE); 
    }
    
    __forceinline ssei_m& operator=( const __m512i& other ) { 
      assert((size_t)this % 16 == 0); 
      _mm512_mask_extpackstorelo_epi32(this,0xf,other,_MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE); 
      return *this; 
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline const int& operator []( const size_t index ) const { assert(index < 4); return i[index]; }
    __forceinline       int& operator []( const size_t index )       { assert(index < 4); return i[index]; }
  };
    
  /*! 4-wide SSE integer type emulated with 16-wide vectors. */
  struct __align(64) ssei_t 
  {
    typedef sseb_t Mask;
    typedef ssei_t Int;
    typedef ssef_t Float;
    enum { size = 4 };  // number of SIMD elements
    union {
      __m512i m512; 
      float  f[4];
      int    i[4];
    };
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline ssei_t           ( ) {}
    __forceinline ssei_t           ( const ssei_t& other ) { m512 = other.m512; }
    __forceinline ssei_t& operator=( const ssei_t& other ) { m512 = other.m512; return *this; }
      
    __forceinline ssei_t           ( const ssei_m& other ) { m512 = _mm512_extload_epi32(&other,_MM_UPCONV_EPI32_NONE,_MM_BROADCAST_4X16,_MM_HINT_NONE); }
    __forceinline ssei_t& operator=( const ssei_m& other ) { m512 = _mm512_extload_epi32(&other,_MM_UPCONV_EPI32_NONE,_MM_BROADCAST_4X16,_MM_HINT_NONE); return *this; }
    
    __forceinline ssei_t& operator=( const int& a ) { m512 = _mm512_extload_epi32(&a,_MM_UPCONV_EPI32_NONE,_MM_BROADCAST_1X16,_MM_HINT_NONE); return *this; }
    
    __forceinline ssei_t( const __m512i a ) : m512(a) {}
    __forceinline operator const __m512i&( void ) const { return m512; }
    __forceinline operator       __m512i&( void )       { return m512; }
    
    __forceinline explicit ssei_t  ( const int* const a ) : m512(_mm512_load_epi32(a)) {}
    __forceinline ssei_t           ( const int& a ) : m512(_mm512_extload_epi32(&a,_MM_UPCONV_EPI32_NONE,_MM_BROADCAST_1X16,_MM_HINT_NONE)) {}
    __forceinline ssei_t           ( int  a, int  b) : m512(ssei_t(ssei_m(a, a, b, b))) {}
    __forceinline ssei_t           ( int  a, int  b, int  c, int  d) : m512(ssei_t(ssei_m(a, b, c, d))) {}
    
    __forceinline explicit ssei_t( const __m512 a ) : m512(_mm512_cvtfxpnt_round_adjustps_epi32(a,_MM_ROUND_MODE_NEAREST,_MM_EXPADJ_NONE)) {}
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline ssei_t( ZeroTy   ) : m512(_mm512_setzero_epi32()) {}
    __forceinline ssei_t( OneTy    ) : m512(_mm512_set1_epi32(1.0f)) {}
    __forceinline ssei_t( PosInfTy ) : m512(_mm512_set1_epi32(pos_inf)) {}
    __forceinline ssei_t( NegInfTy ) : m512(_mm512_set1_epi32(neg_inf)) {}
    __forceinline ssei_t( StepTy   ) : m512(ssei_t(ssei_m(0, 1, 2, 3))) {}
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline const int& operator []( const size_t index ) const { assert(index < 4); return i[index]; }
    __forceinline       int& operator []( const size_t index )       { assert(index < 4); return i[index]; }
  };
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Unary Operators
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline const ssei_t operator +( const ssei_t& a ) { return a; }
  __forceinline const ssei_t operator -( const ssei_t& a ) { return _mm512_sub_epi32(_mm512_setzero_epi32(), a.m512); }
  //__forceinline const ssei_t abs       ( const ssei_t& a ) { return _mm512_absdiff_epi32(a.m512,_mm512_setzero_epi32()); }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline const ssei_t operator +( const ssei_t& a, const ssei_t& b ) { return _mm512_add_epi32(a.m512, b.m512); }
  __forceinline const ssei_t operator +( const ssei_t& a, const int32&  b ) { return a + ssei_t(b); }
  __forceinline const ssei_t operator +( const int32&  a, const ssei_t& b ) { return ssei_t(a) + b; }
  
  __forceinline const ssei_t operator -( const ssei_t& a, const ssei_t& b ) { return _mm512_sub_epi32(a.m512, b.m512); }
  __forceinline const ssei_t operator -( const ssei_t& a, const int32&  b ) { return a - ssei_t(b); }
  __forceinline const ssei_t operator -( const int32&  a, const ssei_t& b ) { return ssei_t(a) - b; }
  
  __forceinline const ssei_t operator *( const ssei_t& a, const ssei_t& b ) { return _mm512_mullo_epi32(a.m512, b.m512); }
  __forceinline const ssei_t operator *( const ssei_t& a, const int32&  b ) { return a * ssei_t(b); }
  __forceinline const ssei_t operator *( const int32&  a, const ssei_t& b ) { return ssei_t(a) * b; }
  
  __forceinline const ssei_t operator &( const ssei_t& a, const ssei_t& b ) { return _mm512_and_epi32(a.m512, b.m512); }
  __forceinline const ssei_t operator &( const ssei_t& a, const int32&  b ) { return a & ssei_t(b); }
  __forceinline const ssei_t operator &( const int32& a, const ssei_t& b ) { return ssei_t(a) & b; }
  
  __forceinline const ssei_t operator |( const ssei_t& a, const ssei_t& b ) { return _mm512_or_epi32(a.m512, b.m512); }
  __forceinline const ssei_t operator |( const ssei_t& a, const int32&  b ) { return a | ssei_t(b); }
  __forceinline const ssei_t operator |( const int32& a, const ssei_t& b ) { return ssei_t(a) | b; }
  
  __forceinline const ssei_t operator ^( const ssei_t& a, const ssei_t& b ) { return _mm512_xor_epi32(a.m512, b.m512); }
  __forceinline const ssei_t operator ^( const ssei_t& a, const int32&  b ) { return a ^ ssei_t(b); }
  __forceinline const ssei_t operator ^( const int32& a, const ssei_t& b ) { return ssei_t(a) ^ b; }
  
  __forceinline const ssei_t operator <<( const ssei_t& a, const int32& n ) { return _mm512_slli_epi32(a.m512, n); }
  __forceinline const ssei_t operator >>( const ssei_t& a, const int32& n ) { return _mm512_srai_epi32(a.m512, n); }
  
  __forceinline const ssei_t sra ( const ssei_t& a, const int32& b ) { return _mm512_srai_epi32(a.m512, b); }
  __forceinline const ssei_t srl ( const ssei_t& a, const int32& b ) { return _mm512_srli_epi32(a.m512, b); }
  
  __forceinline const ssei_t min( const ssei_t& a, const ssei_t& b ) { return _mm512_min_epi32(a.m512, b.m512); }
  __forceinline const ssei_t min( const ssei_t& a, const int32&  b ) { return min(a,ssei_t(b)); }
  __forceinline const ssei_t min( const int32&  a, const ssei_t& b ) { return min(ssei_t(a),b); }
  
  __forceinline const ssei_t max( const ssei_t& a, const ssei_t& b ) { return _mm512_max_epi32(a.m512, b.m512); }
  __forceinline const ssei_t max( const ssei_t& a, const int32&  b ) { return max(a,ssei_t(b)); }
  __forceinline const ssei_t max( const int32&  a, const ssei_t& b ) { return max(ssei_t(a),b); }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline ssei_t& operator +=( ssei_t& a, const ssei_t& b ) { return a = a + b; }
  __forceinline ssei_t& operator +=( ssei_t& a, const int32&  b ) { return a = a + b; }
  
  __forceinline ssei_t& operator -=( ssei_t& a, const ssei_t& b ) { return a = a - b; }
  __forceinline ssei_t& operator -=( ssei_t& a, const int32&  b ) { return a = a - b; }
  
  __forceinline ssei_t& operator *=( ssei_t& a, const ssei_t& b ) { return a = a * b; }
  __forceinline ssei_t& operator *=( ssei_t& a, const int32&  b ) { return a = a * b; }
  
  __forceinline ssei_t& operator &=( ssei_t& a, const ssei_t& b ) { return a = a & b; }
  __forceinline ssei_t& operator &=( ssei_t& a, const int32&  b ) { return a = a & b; }
  
  __forceinline ssei_t& operator |=( ssei_t& a, const ssei_t& b ) { return a = a | b; }
  __forceinline ssei_t& operator |=( ssei_t& a, const int32&  b ) { return a = a | b; }
  
  __forceinline ssei_t& operator <<=( ssei_t& a, const int32&  b ) { return a = a << b; }
  __forceinline ssei_t& operator >>=( ssei_t& a, const int32&  b ) { return a = a >> b; }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators + Select
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline const sseb_t operator ==( const ssei_t& a, const ssei_t& b ) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_EQ); }
  __forceinline const sseb_t operator ==( const ssei_t& a, const int32& b ) { return a == ssei_t(b); }
  __forceinline const sseb_t operator ==( const int32& a, const ssei_t& b ) { return ssei_t(a) == b; }
  
  __forceinline const sseb_t operator !=( const ssei_t& a, const ssei_t& b ) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_NE); }
  __forceinline const sseb_t operator !=( const ssei_t& a, const int32& b ) { return a != ssei_t(b); }
  __forceinline const sseb_t operator !=( const int32& a, const ssei_t& b ) { return ssei_t(a) != b; }
  
  __forceinline const sseb_t operator < ( const ssei_t& a, const ssei_t& b ) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_LT); }
  __forceinline const sseb_t operator < ( const ssei_t& a, const int32& b ) { return a <  ssei_t(b); }
  __forceinline const sseb_t operator < ( const int32& a, const ssei_t& b ) { return ssei_t(a) <  b; }
  
  __forceinline const sseb_t operator >=( const ssei_t& a, const ssei_t& b ) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_GE); }
  __forceinline const sseb_t operator >=( const ssei_t& a, const int32& b ) { return a >= ssei_t(b); }
  __forceinline const sseb_t operator >=( const int32& a, const ssei_t& b ) { return ssei_t(a) >= b; }
  
  __forceinline const sseb_t operator > ( const ssei_t& a, const ssei_t& b ) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_GT); }
  __forceinline const sseb_t operator > ( const ssei_t& a, const int32& b ) { return a >  ssei_t(b); }
  __forceinline const sseb_t operator > ( const int32& a, const ssei_t& b ) { return ssei_t(a) >  b; }
  
  __forceinline const sseb_t operator <=( const ssei_t& a, const ssei_t& b ) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_LE); }
  __forceinline const sseb_t operator <=( const ssei_t& a, const int32& b ) { return a <= ssei_t(b); }
  __forceinline const sseb_t operator <=( const int32& a, const ssei_t& b ) { return ssei_t(a) <= b; }
  
  __forceinline const ssei_t select( const sseb_t& mask, const ssei_t& t, const ssei_t& f ) { return _mm512_mask_blend_epi32(mask, f, t); }
  
  
  ////////////////////////////////////////////////////////////////////////////////
  // Movement/Shifting/Shuffling Functions
  ////////////////////////////////////////////////////////////////////////////////
  
  //__forceinline ssei_t unpacklo( const ssei_t& a, const ssei_t& b ) { return _mm512_castps_si128(_mm512_unpacklo_ps(_mm512_castsi128_ps(a.m512), _mm512_castsi128_ps(b.m512))); }
  //__forceinline ssei_t unpackhi( const ssei_t& a, const ssei_t& b ) { return _mm512_castps_si128(_mm512_unpackhi_ps(_mm512_castsi128_ps(a.m512), _mm512_castsi128_ps(b.m512))); }
  
  template<size_t i0, size_t i1, size_t i2, size_t i3> __forceinline const ssei_t shuffle( const ssei_t& a ) {
    return _mm512_permute4f128_epi32(a, (_MM_PERM_ENUM)_MM_SHUFFLE(i3, i2, i1, i0));
  }
  
  //template<size_t i0, size_t i1, size_t i2, size_t i3> __forceinline const ssei_t shuffle( const ssei_t& a, const ssei_t& b ) {
  //return _mm512_castps_si128(_mm512_shuffle_ps(_mm512_castsi128_ps(a), _mm512_castsi128_ps(b), _MM512_SHUFFLE(i3, i2, i1, i0)));
  //}
  
  template<size_t src> __forceinline int    extract( const ssei_t& a                ) { return a[src]; }
  template<size_t dst> __forceinline ssei_t insert ( const ssei_t& a, const float b ) { ssei_t c = a; c[dst] = b; return c; }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Reductions
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline const ssei_t vreduce_min(const ssei_t& v) { ssei_t h = min(shuffle<1,0,3,2>(v),v); return min(shuffle<2,3,0,1>(h),h); }
  __forceinline const ssei_t vreduce_max(const ssei_t& v) { ssei_t h = max(shuffle<1,0,3,2>(v),v); return max(shuffle<2,3,0,1>(h),h); }
  __forceinline const ssei_t vreduce_add(const ssei_t& v) { ssei_t h = shuffle<1,0,3,2>(v)   + v ; return shuffle<2,3,0,1>(h)   + h ; }
  
  __forceinline int reduce_min(const ssei_t& v) { return extract<0>(vreduce_min(v)); }
  __forceinline int reduce_max(const ssei_t& v) { return extract<0>(vreduce_max(v)); }
  __forceinline int reduce_add(const ssei_t& v) { return extract<0>(vreduce_add(v)); }
  
  __forceinline size_t select_min(const ssei_t& v) { return __bsf(movemask(v == vreduce_min(v))); }
  __forceinline size_t select_max(const ssei_t& v) { return __bsf(movemask(v == vreduce_max(v))); }
  
  __forceinline size_t select_min(const sseb_t& valid, const ssei_t& v) { const ssei_t a = select(valid,v,ssei_t(pos_inf)); return __bsf(movemask(valid & (a == vreduce_min(a)))); }
  __forceinline size_t select_max(const sseb_t& valid, const ssei_t& v) { const ssei_t a = select(valid,v,ssei_t(neg_inf)); return __bsf(movemask(valid & (a == vreduce_max(a)))); }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Output Operators
  ////////////////////////////////////////////////////////////////////////////////
  
  inline std::ostream& operator<<(std::ostream& cout, const ssei_t& a) {
    return cout << "<" << a[0] << ", " << a[1] << ", " << a[2] << ", " << a[3] << ">";
  }
}

#endif


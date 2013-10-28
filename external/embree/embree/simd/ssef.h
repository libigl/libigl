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

#ifndef __EMBREE_SSEF_H__
#define __EMBREE_SSEF_H__

namespace embree
{
  /*! 4-wide SSE float type. */
  struct ssef
  {
    typedef sseb Mask;                    // mask type
    typedef ssei Int;                     // int type
    typedef ssef Float;                   // float type
    
    enum   { size = 4 };  // number of SIMD elements
    union { __m128 m128; float f[4]; int i[4]; }; // data

    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline ssef           ( ) {}
    __forceinline ssef           ( const ssef& other ) { m128 = other.m128; }
    __forceinline ssef& operator=( const ssef& other ) { m128 = other.m128; return *this; }

    __forceinline ssef( const __m128 a ) : m128(a) {}
    __forceinline operator const __m128&( void ) const { return m128; }
    __forceinline operator       __m128&( void )       { return m128; }

    __forceinline explicit ssef( const char* const a ) : m128(_mm_loadu_ps((const float*)a)) {}
    //__forceinline ssef           ( float  a ) : m128(_mm_set1_ps(a)) {}
    __forceinline ssef           ( const float&  a ) : m128(_mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(_mm_load_ss(&a)), _MM_SHUFFLE(0, 0, 0, 0)))) {}
    __forceinline ssef           ( float  a, float  b) : m128(_mm_set_ps(b, a, b, a)) {}
    __forceinline ssef           ( float  a, float  b, float  c, float  d) : m128(_mm_set_ps(d, c, b, a)) {}

    __forceinline explicit ssef( const __m128i a ) : m128(_mm_cvtepi32_ps(a)) {}

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline ssef( ZeroTy   ) : m128(_mm_setzero_ps()) {}
    __forceinline ssef( OneTy    ) : m128(_mm_set_ps(1.0f,1.0f,1.0f,1.0f)) {}
    __forceinline ssef( PosInfTy ) : m128(_mm_set_ps(pos_inf,pos_inf,pos_inf,pos_inf)) {}
    __forceinline ssef( NegInfTy ) : m128(_mm_set_ps(neg_inf,neg_inf,neg_inf,neg_inf)) {}
    __forceinline ssef( StepTy   ) : m128(_mm_set_ps(3.0f, 2.0f, 1.0f, 0.0f)) {}
    __forceinline ssef( NaNTy    ) : m128(_mm_set1_ps(nan)) {}

    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline const float& operator []( const size_t i ) const { assert(i < 4); return f[i]; }
    __forceinline       float& operator []( const size_t i )       { assert(i < 4); return f[i]; }
  };


  ////////////////////////////////////////////////////////////////////////////////
  /// Unary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const ssef operator +( const ssef& a ) { return a; }
  __forceinline const ssef operator -( const ssef& a ) { return _mm_xor_ps(a.m128, _mm_castsi128_ps(_mm_set1_epi32(0x80000000))); }
  __forceinline const ssef abs       ( const ssef& a ) { return _mm_and_ps(a.m128, _mm_castsi128_ps(_mm_set1_epi32(0x7fffffff))); }
  __forceinline const ssef sign      ( const ssef& a ) { return _mm_blendv_ps(ssef(one), -ssef(one), _mm_cmplt_ps (a,ssef(zero))); }
  __forceinline const ssef signmsk   ( const ssef& a ) { return _mm_and_ps(a.m128,_mm_castsi128_ps(_mm_set1_epi32(0x80000000))); }
  
  __forceinline const ssef rcp  ( const ssef& a ) {
    const ssef r = _mm_rcp_ps(a.m128);
    return _mm_sub_ps(_mm_add_ps(r, r), _mm_mul_ps(_mm_mul_ps(r, r), a));
  }
  __forceinline const ssef sqr  ( const ssef& a ) { return _mm_mul_ps(a,a); }
  __forceinline const ssef sqrt ( const ssef& a ) { return _mm_sqrt_ps(a.m128); }
  __forceinline const ssef rsqrt( const ssef& a ) {
    const ssef r = _mm_rsqrt_ps(a.m128);
    return _mm_add_ps(_mm_mul_ps(_mm_set_ps(1.5f, 1.5f, 1.5f, 1.5f), r),
                      _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(a, _mm_set_ps(-0.5f, -0.5f, -0.5f, -0.5f)), r), _mm_mul_ps(r, r)));
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const ssef operator +( const ssef& a, const ssef& b ) { return _mm_add_ps(a.m128, b.m128); }
  __forceinline const ssef operator +( const ssef& a, const float& b ) { return a + ssef(b); }
  __forceinline const ssef operator +( const float& a, const ssef& b ) { return ssef(a) + b; }

  __forceinline const ssef operator -( const ssef& a, const ssef& b ) { return _mm_sub_ps(a.m128, b.m128); }
  __forceinline const ssef operator -( const ssef& a, const float& b ) { return a - ssef(b); }
  __forceinline const ssef operator -( const float& a, const ssef& b ) { return ssef(a) - b; }

  __forceinline const ssef operator *( const ssef& a, const ssef& b ) { return _mm_mul_ps(a.m128, b.m128); }
  __forceinline const ssef operator *( const ssef& a, const float& b ) { return a * ssef(b); }
  __forceinline const ssef operator *( const float& a, const ssef& b ) { return ssef(a) * b; }

  __forceinline const ssef operator /( const ssef& a, const ssef& b ) { return _mm_div_ps(a.m128,b.m128); }
  __forceinline const ssef operator /( const ssef& a, const float& b ) { return a/ssef(b); }
  __forceinline const ssef operator /( const float& a, const ssef& b ) { return ssef(a)/b; }

  __forceinline const ssef operator^( const ssef& a, const ssef& b ) { return _mm_xor_ps(a.m128,b.m128); }
  __forceinline const ssef operator^( const ssef& a, const ssei& b ) { return _mm_xor_ps(a.m128,_mm_castsi128_ps(b.m128)); }

  __forceinline const ssef min( const ssef& a, const ssef& b ) { return _mm_min_ps(a.m128,b.m128); }
  __forceinline const ssef min( const ssef& a, const float& b ) { return _mm_min_ps(a.m128,ssef(b)); }
  __forceinline const ssef min( const float& a, const ssef& b ) { return _mm_min_ps(ssef(a),b.m128); }

  __forceinline const ssef max( const ssef& a, const ssef& b ) { return _mm_max_ps(a.m128,b.m128); }
  __forceinline const ssef max( const ssef& a, const float& b ) { return _mm_max_ps(a.m128,ssef(b)); }
  __forceinline const ssef max( const float& a, const ssef& b ) { return _mm_max_ps(ssef(a),b.m128); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Ternary Operators
  ////////////////////////////////////////////////////////////////////////////////

#if defined(__AVX2__)
  __forceinline const ssef madd  ( const ssef& a, const ssef& b, const ssef& c) { return _mm_fmadd_ps(a,b,c); }
  __forceinline const ssef msub  ( const ssef& a, const ssef& b, const ssef& c) { return _mm_fmsub_ps(a,b,c); }
  __forceinline const ssef nmadd ( const ssef& a, const ssef& b, const ssef& c) { return _mm_fnmadd_ps(a,b,c); }
  __forceinline const ssef nmsub ( const ssef& a, const ssef& b, const ssef& c) { return _mm_fnmsub_ps(a,b,c); }
#else
  __forceinline const ssef madd  ( const ssef& a, const ssef& b, const ssef& c) { return a*b+c; }
  __forceinline const ssef msub  ( const ssef& a, const ssef& b, const ssef& c) { return a*b-c; }
  __forceinline const ssef nmadd ( const ssef& a, const ssef& b, const ssef& c) { return -a*b-c;}
  __forceinline const ssef nmsub ( const ssef& a, const ssef& b, const ssef& c) { return c-a*b; }
#endif

  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline ssef& operator +=( ssef& a, const ssef& b ) { return a = a + b; }
  __forceinline ssef& operator +=( ssef& a, const float& b ) { return a = a + b; }

  __forceinline ssef& operator -=( ssef& a, const ssef& b ) { return a = a - b; }
  __forceinline ssef& operator -=( ssef& a, const float& b ) { return a = a - b; }

  __forceinline ssef& operator *=( ssef& a, const ssef& b ) { return a = a * b; }
  __forceinline ssef& operator *=( ssef& a, const float& b ) { return a = a * b; }

  __forceinline ssef& operator /=( ssef& a, const ssef& b ) { return a = a / b; }
  __forceinline ssef& operator /=( ssef& a, const float& b ) { return a = a / b; }

  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators + Select
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const sseb operator ==( const ssef& a, const ssef& b ) { return _mm_cmpeq_ps (a.m128, b.m128); }
  __forceinline const sseb operator ==( const ssef& a, const float& b ) { return a == ssef(b); }
  __forceinline const sseb operator ==( const float& a, const ssef& b ) { return ssef(a) == b; }

  __forceinline const sseb operator !=( const ssef& a, const ssef& b ) { return _mm_cmpneq_ps(a.m128, b.m128); }
  __forceinline const sseb operator !=( const ssef& a, const float& b ) { return a != ssef(b); }
  __forceinline const sseb operator !=( const float& a, const ssef& b ) { return ssef(a) != b; }

  __forceinline const sseb operator < ( const ssef& a, const ssef& b ) { return _mm_cmplt_ps (a.m128, b.m128); }
  __forceinline const sseb operator < ( const ssef& a, const float& b ) { return a <  ssef(b); }
  __forceinline const sseb operator < ( const float& a, const ssef& b ) { return ssef(a) <  b; }

  __forceinline const sseb operator >=( const ssef& a, const ssef& b ) { return _mm_cmpnlt_ps(a.m128, b.m128); }
  __forceinline const sseb operator >=( const ssef& a, const float& b ) { return a >= ssef(b); }
  __forceinline const sseb operator >=( const float& a, const ssef& b ) { return ssef(a) >= b; }

  __forceinline const sseb operator > ( const ssef& a, const ssef& b ) { return _mm_cmpnle_ps(a.m128, b.m128); }
  __forceinline const sseb operator > ( const ssef& a, const float& b ) { return a >  ssef(b); }
  __forceinline const sseb operator > ( const float& a, const ssef& b ) { return ssef(a) >  b; }

  __forceinline const sseb operator <=( const ssef& a, const ssef& b ) { return _mm_cmple_ps (a.m128, b.m128); }
  __forceinline const sseb operator <=( const ssef& a, const float& b ) { return a <= ssef(b); }
  __forceinline const sseb operator <=( const float& a, const ssef& b ) { return ssef(a) <= b; }

 __forceinline const ssef select( const sseb& mask, const ssef& t, const ssef& f ) { 
    return _mm_blendv_ps(f, t, mask); 
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Rounding Functions
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const ssef round_even( const ssef& a ) { return _mm_round_ps(a, _MM_FROUND_TO_NEAREST_INT); }
  __forceinline const ssef round_down( const ssef& a ) { return _mm_round_ps(a, _MM_FROUND_TO_NEG_INF    ); }
  __forceinline const ssef round_up  ( const ssef& a ) { return _mm_round_ps(a, _MM_FROUND_TO_POS_INF    ); }
  __forceinline const ssef round_zero( const ssef& a ) { return _mm_round_ps(a, _MM_FROUND_TO_ZERO       ); }
  __forceinline const ssef floor     ( const ssef& a ) { return _mm_round_ps(a, _MM_FROUND_TO_NEG_INF    ); }
  __forceinline const ssef ceil      ( const ssef& a ) { return _mm_round_ps(a, _MM_FROUND_TO_POS_INF    ); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Movement/Shifting/Shuffling Functions
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline ssef unpacklo( const ssef& a, const ssef& b ) { return _mm_unpacklo_ps(a.m128, b.m128); }
  __forceinline ssef unpackhi( const ssef& a, const ssef& b ) { return _mm_unpackhi_ps(a.m128, b.m128); }

  template<size_t i0, size_t i1, size_t i2, size_t i3> __forceinline const ssef shuffle( const ssef& b ) {
    return _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(b), _MM_SHUFFLE(i3, i2, i1, i0)));
  }

  template<size_t i0, size_t i1, size_t i2, size_t i3> __forceinline const ssef shuffle( const ssef& a, const ssef& b ) {
    return _mm_shuffle_ps(a, b, _MM_SHUFFLE(i3, i2, i1, i0));
  }

  __forceinline const ssef shuffle8(const ssef& a, const ssei& shuf) { 
    return _mm_castsi128_ps(_mm_shuffle_epi8(_mm_castps_si128(a), shuf)); 
  }

  template<> __forceinline const ssef shuffle<0, 0, 2, 2>( const ssef& b ) { return _mm_moveldup_ps(b); }
  template<> __forceinline const ssef shuffle<1, 1, 3, 3>( const ssef& b ) { return _mm_movehdup_ps(b); }
  template<> __forceinline const ssef shuffle<0, 1, 0, 1>( const ssef& b ) { return _mm_castpd_ps(_mm_movedup_pd(_mm_castps_pd(b))); }

  template<size_t i> __forceinline float extract( const ssef& a ) { return _mm_cvtss_f32(shuffle<i,i,i,i>(a)); }
  template<size_t dst, size_t src, size_t clr> __forceinline const ssef insert( const ssef& a, const ssef& b ) { return _mm_insert_ps(a, b, (dst << 4) | (src << 6) | clr); }
  template<size_t dst, size_t src> __forceinline const ssef insert( const ssef& a, const ssef& b ) { return insert<dst, src, 0>(a, b); }
  template<size_t dst>             __forceinline const ssef insert( const ssef& a, const float b ) { return insert<dst,      0>(a, _mm_set_ss(b)); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Transpose
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline void transpose(const ssef& r0, const ssef& r1, const ssef& r2, const ssef& r3, ssef& c0, ssef& c1, ssef& c2, ssef& c3) {
    ssef l02 = unpacklo(r0,r2);
    ssef h02 = unpackhi(r0,r2);
    ssef l13 = unpacklo(r1,r3);
    ssef h13 = unpackhi(r1,r3);
    c0 = unpacklo(l02,l13);
    c1 = unpackhi(l02,l13);
    c2 = unpacklo(h02,h13);
    c3 = unpackhi(h02,h13);
  }

  __forceinline void transpose(const ssef& r0, const ssef& r1, const ssef& r2, const ssef& r3, ssef& c0, ssef& c1, ssef& c2) {
    ssef l02 = unpacklo(r0,r2);
    ssef h02 = unpackhi(r0,r2);
    ssef l13 = unpacklo(r1,r3);
    ssef h13 = unpackhi(r1,r3);
    c0 = unpacklo(l02,l13);
    c1 = unpackhi(l02,l13);
    c2 = unpacklo(h02,h13);
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Reductions
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const ssef vreduce_min(const ssef& v) { ssef h = min(shuffle<1,0,3,2>(v),v); return min(shuffle<2,3,0,1>(h),h); }
  __forceinline const ssef vreduce_max(const ssef& v) { ssef h = max(shuffle<1,0,3,2>(v),v); return max(shuffle<2,3,0,1>(h),h); }
  __forceinline const ssef vreduce_add(const ssef& v) { ssef h = shuffle<1,0,3,2>(v)   + v ; return shuffle<2,3,0,1>(h)   + h ; }

  __forceinline float reduce_min(const ssef& v) { return _mm_cvtss_f32(vreduce_min(v)); }
  __forceinline float reduce_max(const ssef& v) { return _mm_cvtss_f32(vreduce_max(v)); }
  __forceinline float reduce_add(const ssef& v) { return _mm_cvtss_f32(vreduce_add(v)); }

  __forceinline size_t select_min(const ssef& v) { return __bsf(movemask(v == vreduce_min(v))); }
  __forceinline size_t select_max(const ssef& v) { return __bsf(movemask(v == vreduce_max(v))); }

  __forceinline size_t select_min(const sseb& valid, const ssef& v) { const ssef a = select(valid,v,ssef(pos_inf)); return __bsf(movemask(valid & (a == vreduce_min(a)))); }
  __forceinline size_t select_max(const sseb& valid, const ssef& v) { const ssef a = select(valid,v,ssef(neg_inf)); return __bsf(movemask(valid & (a == vreduce_max(a)))); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Output Operators
  ////////////////////////////////////////////////////////////////////////////////

  inline std::ostream& operator<<(std::ostream& cout, const ssef& a) {
    return cout << "<" << a[0] << ", " << a[1] << ", " << a[2] << ", " << a[3] << ">";
  }

  typedef ssef ssef_m;
  typedef ssef ssef_t;
}

#endif

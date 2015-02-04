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

#include "simd/sse.h"
#include "math.h"

namespace embree
{
  ////////////////////////////////////////////////////////////////////////////////
  /// SSE Vec3fa Type
  ////////////////////////////////////////////////////////////////////////////////

  struct __aligned(16) Vec3fa
  {
    typedef float Scalar;
    enum { N = 3 };
    union {
      __m128 m128;
      struct { float x,y,z; union { int a; unsigned u; float w; }; };
    };

    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vec3fa( ) {}
    __forceinline Vec3fa( const __m128 a ) : m128(a) {}

    __forceinline Vec3fa            ( const Vec3<float>& other  ) { x = other.x; y = other.y; z = other.z; }
    __forceinline Vec3fa& operator =( const Vec3<float>& other ) { x = other.x; y = other.y; z = other.z; return *this; }

    __forceinline Vec3fa            ( const Vec3fa& other ) { m128 = other.m128; }
    __forceinline Vec3fa& operator =( const Vec3fa& other ) { m128 = other.m128; return *this; }

    __forceinline Vec3fa( const float a ) : m128(_mm_set1_ps(a)) {}
    __forceinline Vec3fa( const float x, const float y, const float z) : m128(_mm_set_ps(z, z, y, x)) {}

    __forceinline explicit Vec3fa( const Vec3fa& other, const int      a1) { m128 = other.m128; a = a1; }
    __forceinline explicit Vec3fa( const Vec3fa& other, const unsigned a1) { m128 = other.m128; a = a1; }
    __forceinline explicit Vec3fa( const Vec3fa& other, const float    w1) { m128 = other.m128; w = w1; }
    __forceinline explicit Vec3fa( const float x, const float y, const float z, const int      a) : x(x), y(y), z(z), a(a) {}
    __forceinline explicit Vec3fa( const float x, const float y, const float z, const unsigned a) : x(x), y(y), z(z), a(a) {}
    __forceinline explicit Vec3fa( const float x, const float y, const float z, const float    w) : x(x), y(y), z(z), w(w) {}

    __forceinline explicit Vec3fa( const __m128i a ) : m128(_mm_cvtepi32_ps(a)) {}

    __forceinline operator const __m128&( void ) const { return m128; }
    __forceinline operator       __m128&( void )       { return m128; }

#if defined (__SSE4_1__)
    friend __forceinline const Vec3fa copy_a( const Vec3fa& a, const Vec3fa& b ) { return _mm_insert_ps(a, b, (3 << 4) | (3 << 6)); }
#else
    friend __forceinline const Vec3fa copy_a( const Vec3fa& a, const Vec3fa& b ) { Vec3fa c = a; c.a = b.a; return c; }
#endif

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vec3fa( ZeroTy   ) : m128(_mm_setzero_ps()) {}
    __forceinline Vec3fa( OneTy    ) : m128(_mm_set1_ps(1.0f)) {}
    __forceinline Vec3fa( PosInfTy ) : m128(_mm_set1_ps(pos_inf)) {}
    __forceinline Vec3fa( NegInfTy ) : m128(_mm_set1_ps(neg_inf)) {}

    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline const float& operator []( const size_t index ) const { assert(index < 3); return (&x)[index]; }
    __forceinline       float& operator []( const size_t index )       { assert(index < 3); return (&x)[index]; }
  };
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Unary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vec3fa operator +( const Vec3fa& a ) { return a; }
  __forceinline const Vec3fa operator -( const Vec3fa& a ) {
    const __m128 mask = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
    return _mm_xor_ps(a.m128, mask);
  }
  __forceinline const Vec3fa abs  ( const Vec3fa& a ) {
    const __m128 mask = _mm_castsi128_ps(_mm_set1_epi32(0x7fffffff));
    return _mm_and_ps(a.m128, mask);
  }
  __forceinline const Vec3fa sign ( const Vec3fa& a ) {
    return blendv_ps(Vec3fa(one), -Vec3fa(one), _mm_cmplt_ps (a,Vec3fa(zero))); 
  }
  __forceinline const Vec3fa rcp  ( const Vec3fa& a ) {
    const Vec3fa r = _mm_rcp_ps(a.m128);
    return _mm_sub_ps(_mm_add_ps(r, r), _mm_mul_ps(_mm_mul_ps(r, r), a));
  }
  __forceinline const Vec3fa sqrt ( const Vec3fa& a ) { return _mm_sqrt_ps(a.m128); }
  __forceinline const Vec3fa sqr  ( const Vec3fa& a ) { return _mm_mul_ps(a,a); }
  __forceinline const Vec3fa rsqrt( const Vec3fa& a ) {
    __m128 r = _mm_rsqrt_ps(a.m128);
    return _mm_add_ps(_mm_mul_ps(_mm_set1_ps(1.5f),r), _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(a, _mm_set1_ps(-0.5f)), r), _mm_mul_ps(r, r)));
  }
  __forceinline const Vec3fa zero_fix(const Vec3fa& a) { return blendv_ps(a, _mm_set1_ps(1E-10f), _mm_cmpeq_ps (a.m128, _mm_setzero_ps())); }
  __forceinline const Vec3fa rcp_safe(const Vec3fa& a) { return rcp(zero_fix(a)); }

  __forceinline Vec3fa log ( const Vec3fa& a ) { 
    return Vec3fa(logf(a.x),logf(a.y),logf(a.z));
  }

  __forceinline Vec3fa exp ( const Vec3fa& a ) { 
    return Vec3fa(expf(a.x),expf(a.y),expf(a.z));
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vec3fa operator +( const Vec3fa& a, const Vec3fa& b ) { return _mm_add_ps(a.m128, b.m128); }
  __forceinline const Vec3fa operator -( const Vec3fa& a, const Vec3fa& b ) { return _mm_sub_ps(a.m128, b.m128); }
  __forceinline const Vec3fa operator *( const Vec3fa& a, const Vec3fa& b ) { return _mm_mul_ps(a.m128, b.m128); }
  __forceinline const Vec3fa operator *( const Vec3fa& a, const float b ) { return a * Vec3fa(b); }
  __forceinline const Vec3fa operator *( const float a, const Vec3fa& b ) { return Vec3fa(a) * b; }
  __forceinline const Vec3fa operator /( const Vec3fa& a, const Vec3fa& b ) { return _mm_div_ps(a.m128,b.m128); }
  __forceinline const Vec3fa operator /( const Vec3fa& a, const float b        ) { return _mm_div_ps(a.m128,_mm_set1_ps(b)); }
  __forceinline const Vec3fa operator /( const        float a, const Vec3fa& b ) { return _mm_div_ps(_mm_set1_ps(a),b.m128); }

  __forceinline const Vec3fa min( const Vec3fa& a, const Vec3fa& b ) { return _mm_min_ps(a.m128,b.m128); }
  __forceinline const Vec3fa max( const Vec3fa& a, const Vec3fa& b ) { return _mm_max_ps(a.m128,b.m128); }

#if defined(__SSE4_1__)
    __forceinline Vec3fa mini(const Vec3fa& a, const Vec3fa& b) {
      const ssei ai = _mm_castps_si128(a);
      const ssei bi = _mm_castps_si128(b);
      const ssei ci = _mm_min_epi32(ai,bi);
      return _mm_castsi128_ps(ci);
    }
#endif
    
#if defined(__SSE4_1__)
    __forceinline Vec3fa maxi(const Vec3fa& a, const Vec3fa& b) {
      const ssei ai = _mm_castps_si128(a);
      const ssei bi = _mm_castps_si128(b);
      const ssei ci = _mm_max_epi32(ai,bi);
      return _mm_castsi128_ps(ci);
    }
#endif

    __forceinline Vec3fa pow ( const Vec3fa& a, const float& b ) { 
      return Vec3fa(powf(a.x,b),powf(a.y,b),powf(a.z,b));
    }

  ////////////////////////////////////////////////////////////////////////////////
  /// Ternary Operators
  ////////////////////////////////////////////////////////////////////////////////

#if defined(__AVX2__)
  __forceinline Vec3fa madd  ( const Vec3fa& a, const Vec3fa& b, const Vec3fa& c) { return _mm_fmadd_ps(a,b,c); }
  __forceinline Vec3fa msub  ( const Vec3fa& a, const Vec3fa& b, const Vec3fa& c) { return _mm_fmsub_ps(a,b,c); }
  __forceinline Vec3fa nmadd ( const Vec3fa& a, const Vec3fa& b, const Vec3fa& c) { return _mm_fnmadd_ps(a,b,c); }
  __forceinline Vec3fa nmsub ( const Vec3fa& a, const Vec3fa& b, const Vec3fa& c) { return _mm_fnmsub_ps(a,b,c); }
#else
  __forceinline Vec3fa madd  ( const Vec3fa& a, const Vec3fa& b, const Vec3fa& c) { return a*b+c; }
  __forceinline Vec3fa msub  ( const Vec3fa& a, const Vec3fa& b, const Vec3fa& c) { return a*b-c; }
  __forceinline Vec3fa nmadd ( const Vec3fa& a, const Vec3fa& b, const Vec3fa& c) { return -a*b-c;}
  __forceinline Vec3fa nmsub ( const Vec3fa& a, const Vec3fa& b, const Vec3fa& c) { return c-a*b; }
#endif

  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline Vec3fa& operator +=( Vec3fa& a, const Vec3fa& b ) { return a = a + b; }
  __forceinline Vec3fa& operator -=( Vec3fa& a, const Vec3fa& b ) { return a = a - b; }
  __forceinline Vec3fa& operator *=( Vec3fa& a, const Vec3fa& b ) { return a = a * b; }
  __forceinline Vec3fa& operator *=( Vec3fa& a, const float   b ) { return a = a * b; }
  __forceinline Vec3fa& operator /=( Vec3fa& a, const Vec3fa& b ) { return a = a / b; }
  __forceinline Vec3fa& operator /=( Vec3fa& a, const float   b ) { return a = a / b; }

  ////////////////////////////////////////////////////////////////////////////////
  /// Reductions
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline float reduce_add(const Vec3fa& v) { return v.x+v.y+v.z; }
  __forceinline float reduce_mul(const Vec3fa& v) { return v.x*v.y*v.z; }
  __forceinline float reduce_min(const Vec3fa& v) { return min(v.x,v.y,v.z); }
  __forceinline float reduce_max(const Vec3fa& v) { return max(v.x,v.y,v.z); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline bool operator ==( const Vec3fa& a, const Vec3fa& b ) { return (_mm_movemask_ps(_mm_cmpeq_ps (a.m128, b.m128)) & 7) == 7; }
  __forceinline bool operator !=( const Vec3fa& a, const Vec3fa& b ) { return (_mm_movemask_ps(_mm_cmpneq_ps(a.m128, b.m128)) & 7) != 0; }

  __forceinline Vec3ba eq_mask( const Vec3fa& a, const Vec3fa& b ) { return _mm_cmpeq_ps (a.m128, b.m128); }
  __forceinline Vec3ba neq_mask(const Vec3fa& a, const Vec3fa& b ) { return _mm_cmpneq_ps(a.m128, b.m128); }
  __forceinline Vec3ba lt_mask( const Vec3fa& a, const Vec3fa& b ) { return _mm_cmplt_ps (a.m128, b.m128); }
  __forceinline Vec3ba le_mask( const Vec3fa& a, const Vec3fa& b ) { return _mm_cmple_ps (a.m128, b.m128); }
  __forceinline Vec3ba gt_mask( const Vec3fa& a, const Vec3fa& b ) { return _mm_cmpnle_ps(a.m128, b.m128); }
  __forceinline Vec3ba ge_mask( const Vec3fa& a, const Vec3fa& b ) { return _mm_cmpnlt_ps(a.m128, b.m128); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Euclidian Space Operators
  ////////////////////////////////////////////////////////////////////////////////

#if defined(__SSE4_1__)
  __forceinline float dot ( const Vec3fa& a, const Vec3fa& b ) {
    return _mm_cvtss_f32(_mm_dp_ps(a,b,0x7F));
  }
#else
  __forceinline float dot ( const Vec3fa& a, const Vec3fa& b ) {
    return reduce_add(a*b);
  }
#endif

  __forceinline Vec3fa cross ( const Vec3fa& a, const Vec3fa& b ) 
  {
    ssef a0 = ssef(a);
    ssef b0 = shuffle<1,2,0,3>(ssef(b));
    ssef a1 = shuffle<1,2,0,3>(ssef(a));
    ssef b1 = ssef(b);
    return Vec3fa(shuffle<1,2,0,3>(msub(a0,b0,a1*b1)));
  }

  __forceinline float  length   ( const Vec3fa& a )                  { return sqrt(dot(a,a)); }
  __forceinline Vec3fa normalize( const Vec3fa& a )                  { return a*rsqrt(dot(a,a)); }

  __forceinline Vec3fa normalize_safe( const Vec3fa& a ) { 
    const float d = dot(a,a);
    if (unlikely(d == 0.0f)) 
      return a;
    else
      return a*rsqrt(d);
  }

  __forceinline float  distance ( const Vec3fa& a, const Vec3fa& b ) { return length(a-b); }
  __forceinline float  halfArea ( const Vec3fa& d )                  { return d.x*(d.y+d.z)+d.y*d.z; }
  __forceinline Vec3fa reflect  (const Vec3fa& V, const Vec3fa& N)   { return 2.0f*dot(V,N)*N-V; } 

  ////////////////////////////////////////////////////////////////////////////////
  /// Select
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vec3fa select( bool s, const Vec3fa& t, const Vec3fa& f ) {
    __m128 mask = s ? _mm_castsi128_ps(_mm_cmpeq_epi32(_mm_setzero_si128(), _mm_setzero_si128())) : _mm_setzero_ps();
    return blendv_ps(f, t, mask);
  }

  __forceinline const Vec3fa select( const Vec3ba& s, const Vec3fa& t, const Vec3fa& f ) {
    return blendv_ps(f, t, s);
  }

  __forceinline int maxDim ( const Vec3fa& a ) 
  { 
    if (a.x > a.y) {
      if (a.x > a.z) return 0; else return 2;
    } else {
      if (a.y > a.z) return 1; else return 2;
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Rounding Functions
  ////////////////////////////////////////////////////////////////////////////////

#if defined (__SSE4_1__)
  __forceinline const Vec3fa trunc( const Vec3fa& a ) { return _mm_round_ps(a, _MM_FROUND_TO_NEAREST_INT); }
  __forceinline const Vec3fa floor( const Vec3fa& a ) { return _mm_round_ps(a, _MM_FROUND_TO_NEG_INF    ); }
  __forceinline const Vec3fa ceil ( const Vec3fa& a ) { return _mm_round_ps(a, _MM_FROUND_TO_POS_INF    ); }
#else
  __forceinline const Vec3fa trunc( const Vec3fa& a ) { return Vec3fa(truncf(a.x),truncf(a.y),truncf(a.z)); }
  __forceinline const Vec3fa floor( const Vec3fa& a ) { return Vec3fa(floorf(a.x),floorf(a.y),floorf(a.z)); }
  __forceinline const Vec3fa ceil ( const Vec3fa& a ) { return Vec3fa(ceilf (a.x),ceilf (a.y),ceilf (a.z)); }
#endif

  ////////////////////////////////////////////////////////////////////////////////
  /// Output Operators
  ////////////////////////////////////////////////////////////////////////////////

  inline std::ostream& operator<<(std::ostream& cout, const Vec3fa& a) {
    return cout << "(" << a.x << ", " << a.y << ", " << a.z << ")";
  }
  
  typedef Vec3fa Vec3fa_t;
}

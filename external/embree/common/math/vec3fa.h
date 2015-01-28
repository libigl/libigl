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

#ifndef __EMBREE_VECTOR3F_SSE_H__
#define __EMBREE_VECTOR3F_SSE_H__

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
      struct { float x,y,z; union { int a; float w; }; };
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

    __forceinline explicit Vec3fa( const Vec3fa& other, const int a1) { m128 = other.m128; a = a1; }
    __forceinline explicit Vec3fa( const float x, const float y, const float z, const int a) : x(x), y(y), z(z), a(a) {}

    __forceinline explicit Vec3fa( const __m128i a ) : m128(_mm_cvtepi32_ps(a)) {}

    __forceinline operator const __m128&( void ) const { return m128; }
    __forceinline operator       __m128&( void )       { return m128; }

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
    return _mm_blendv_ps(Vec3fa(one), -Vec3fa(one), _mm_cmplt_ps (a,Vec3fa(zero))); 
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
  __forceinline const Vec3fa zero_fix(const Vec3fa& a) { return _mm_blendv_ps(a, _mm_set1_ps(1E-10f), _mm_cmpeq_ps (a.m128, _mm_setzero_ps())); }
  __forceinline const Vec3fa rcp_safe(const Vec3fa& a) { return rcp(zero_fix(a)); }

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

  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline Vec3fa& operator +=( Vec3fa& a, const Vec3fa& b ) { return a = a + b; }
  __forceinline Vec3fa& operator -=( Vec3fa& a, const Vec3fa& b ) { return a = a - b; }
  __forceinline Vec3fa& operator *=( Vec3fa& a, const float b ) { return a = a * b; }
  __forceinline Vec3fa& operator /=( Vec3fa& a, const float b ) { return a = a / b; }

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

  __forceinline float  length   ( const Vec3fa& a )                    { return sqrt(dot(a,a)); }
  __forceinline Vec3fa normalize( const Vec3fa& a )                    { return a*rsqrt(dot(a,a)); }
  __forceinline float  distance ( const Vec3fa& a, const Vec3fa& b ) { return length(a-b); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Select
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vec3fa select( bool s, const Vec3fa& t, const Vec3fa& f ) {
    __m128 mask = s ? _mm_castsi128_ps(_mm_cmpeq_epi32(_mm_setzero_si128(), _mm_setzero_si128())) : _mm_setzero_ps();
    return _mm_blendv_ps(f, t, mask);
  }

  __forceinline const Vec3fa select( const Vec3ba& s, const Vec3fa& t, const Vec3fa& f ) {
    return _mm_blendv_ps(f, t, s);
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
  /// Output Operators
  ////////////////////////////////////////////////////////////////////////////////

  inline std::ostream& operator<<(std::ostream& cout, const Vec3fa& a) {
    return cout << "(" << a.x << ", " << a.y << ", " << a.z << ")";
  }
}

#endif

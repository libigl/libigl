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

#ifndef __EMBREE_VEC3F_SSE_H__
#define __EMBREE_VEC3F_SSE_H__

#include "simd/sse.h"
#include "math.h"

namespace embree
{
  ////////////////////////////////////////////////////////////////////////////////
  /// SSE Vector3f Type
  ////////////////////////////////////////////////////////////////////////////////

  struct Vector3f;

  /* 3 aligned floats as memory representation */
  struct __align(16) Vec3fa
  {
    typedef float Scalar;
    enum { N = 3 };
    union {
      __m128 m128;
      struct { float x,y,z; int a; };
    };
    
    __forceinline Vec3fa ( ) {}
    __forceinline Vec3fa ( float x                   ) : x(x), y(x), z(x), a(0) {}
    __forceinline Vec3fa ( float x, float y, float z ) : x(x), y(y), z(z), a(0) {}

    __forceinline Vec3fa           ( const Vector3f& other );
    __forceinline Vec3fa& operator=( const Vector3f& other );

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

  struct __align(16) Vector3f
  {
    typedef float Scalar;
    enum { N = 3 };
    union {
      __m128 m128;
      struct { float x,y,z; int a; };
    };

    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vector3f( ) {}
    __forceinline Vector3f( const __m128 a ) : m128(a) {}

    __forceinline Vector3f            ( const Vector3f& other ) { m128 = other.m128; }
    __forceinline Vector3f& operator =( const Vector3f& other ) { m128 = other.m128; return *this; }

    __forceinline Vector3f            ( const Vec3fa& other ) { m128 = other.m128; }
    __forceinline Vector3f& operator =( const Vec3fa& other ) { m128 = other.m128; return *this; }

    __forceinline explicit Vector3f( const float a ) : m128(_mm_set1_ps(a)) {}
    __forceinline explicit Vector3f( const float x, const float y, const float z) : m128(_mm_set_ps(z, z, y, x)) {}
    __forceinline explicit Vector3f( const float* const a, const ssize_t stride = 1 ) : m128(_mm_set_ps(0.0f,a[2*stride],a[stride],a[0])) {}
    __forceinline explicit Vector3f( const __m128i a ) : m128(_mm_cvtepi32_ps(a)) {}

    __forceinline operator const __m128&( void ) const { return m128; }
    __forceinline operator       __m128&( void )       { return m128; }

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vector3f( ZeroTy   ) : m128(_mm_setzero_ps()) {}
    __forceinline Vector3f( OneTy    ) : m128(_mm_set1_ps(1.0f)) {}
    __forceinline Vector3f( PosInfTy ) : m128(_mm_set1_ps(pos_inf)) {}
    __forceinline Vector3f( NegInfTy ) : m128(_mm_set1_ps(neg_inf)) {}

    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline const float& operator []( const size_t index ) const { assert(index < 3); return (&x)[index]; }
    __forceinline       float& operator []( const size_t index )       { assert(index < 3); return (&x)[index]; }
  };

  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline          Vec3fa::Vec3fa  ( const Vector3f& other ) { m128 = other.m128; }
  __forceinline Vec3fa& Vec3fa::operator=( const Vector3f& other ) { m128 = other.m128; return *this; }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Unary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vector3f operator +( const Vector3f& a ) { return a; }
  __forceinline const Vector3f operator -( const Vector3f& a ) {
    const __m128 mask = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
    return _mm_xor_ps(a.m128, mask);
  }
  __forceinline const Vector3f abs  ( const Vector3f& a ) {
    const __m128 mask = _mm_castsi128_ps(_mm_set1_epi32(0x7fffffff));
    return _mm_and_ps(a.m128, mask);
  }
  __forceinline const Vector3f sign ( const Vector3f& a ) {
    return _mm_blendv_ps(Vector3f(one), -Vector3f(one), _mm_cmplt_ps (a,Vector3f(zero))); 
  }
  __forceinline const Vector3f rcp  ( const Vector3f& a ) {
    const Vector3f r = _mm_rcp_ps(a.m128);
    return _mm_sub_ps(_mm_add_ps(r, r), _mm_mul_ps(_mm_mul_ps(r, r), a));
  }
  __forceinline const Vector3f sqrt ( const Vector3f& a ) { return _mm_sqrt_ps(a.m128); }
  __forceinline const Vector3f sqr  ( const Vector3f& a ) { return _mm_mul_ps(a,a); }
  __forceinline const Vector3f rsqrt( const Vector3f& a ) {
    __m128 r = _mm_rsqrt_ps(a.m128);
    return _mm_add_ps(_mm_mul_ps(_mm_set1_ps(1.5f),r), _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(a, _mm_set1_ps(-0.5f)), r), _mm_mul_ps(r, r)));
  }
  __forceinline const Vector3f zero_fix(const Vector3f& a) { return _mm_blendv_ps(a, _mm_set1_ps(1E-10f), _mm_cmpeq_ps (a.m128, _mm_setzero_ps())); }
  __forceinline const Vector3f rcp_safe(const Vector3f& a) { return rcp(zero_fix(a)); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vector3f operator +( const Vector3f& a, const Vector3f& b ) { return _mm_add_ps(a.m128, b.m128); }
  __forceinline const Vector3f operator -( const Vector3f& a, const Vector3f& b ) { return _mm_sub_ps(a.m128, b.m128); }
  __forceinline const Vector3f operator *( const Vector3f& a, const Vector3f& b ) { return _mm_mul_ps(a.m128, b.m128); }
  __forceinline const Vector3f operator *( const Vector3f& a, const float b ) { return a * Vector3f(b); }
  __forceinline const Vector3f operator *( const float a, const Vector3f& b ) { return Vector3f(a) * b; }
  __forceinline const Vector3f operator /( const Vector3f& a, const Vector3f& b ) { return _mm_div_ps(a.m128,b.m128); }
  __forceinline const Vector3f operator /( const Vector3f& a, const float b        ) { return _mm_div_ps(a.m128,_mm_set1_ps(b)); }
  __forceinline const Vector3f operator /( const        float a, const Vector3f& b ) { return _mm_div_ps(_mm_set1_ps(a),b.m128); }

  __forceinline const Vector3f min( const Vector3f& a, const Vector3f& b ) { return _mm_min_ps(a.m128,b.m128); }
  __forceinline const Vector3f max( const Vector3f& a, const Vector3f& b ) { return _mm_max_ps(a.m128,b.m128); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline Vector3f& operator +=( Vector3f& a, const Vector3f& b ) { return a = a + b; }
  __forceinline Vector3f& operator -=( Vector3f& a, const Vector3f& b ) { return a = a - b; }
  __forceinline Vector3f& operator *=( Vector3f& a, const float b ) { return a = a * b; }
  __forceinline Vector3f& operator /=( Vector3f& a, const float b ) { return a = a / b; }

  ////////////////////////////////////////////////////////////////////////////////
  /// Reductions
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline float reduce_add(const Vector3f& v) { return v.x+v.y+v.z; }
  __forceinline float reduce_mul(const Vector3f& v) { return v.x*v.y*v.z; }
  __forceinline float reduce_min(const Vector3f& v) { return min(v.x,v.y,v.z); }
  __forceinline float reduce_max(const Vector3f& v) { return max(v.x,v.y,v.z); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline bool operator ==( const Vector3f& a, const Vector3f& b ) { return (_mm_movemask_ps(_mm_cmpeq_ps (a.m128, b.m128)) & 7) == 7; }
  __forceinline bool operator !=( const Vector3f& a, const Vector3f& b ) { return (_mm_movemask_ps(_mm_cmpneq_ps(a.m128, b.m128)) & 7) != 0; }
  __forceinline bool operator < ( const Vector3f& a, const Vector3f& b ) {
    if (a.x != b.x) return a.x < b.x;
    if (a.y != b.y) return a.y < b.y;
    if (a.z != b.z) return a.z < b.z;
    return false;
  }

  __forceinline Vector3b eq_mask( const Vector3f& a, const Vector3f& b ) { return _mm_cmpeq_ps (a.m128, b.m128); }
  __forceinline Vector3b neq_mask(const Vector3f& a, const Vector3f& b ) { return _mm_cmpneq_ps(a.m128, b.m128); }
  __forceinline Vector3b lt_mask( const Vector3f& a, const Vector3f& b ) { return _mm_cmplt_ps (a.m128, b.m128); }
  __forceinline Vector3b le_mask( const Vector3f& a, const Vector3f& b ) { return _mm_cmple_ps (a.m128, b.m128); }
  __forceinline Vector3b gt_mask( const Vector3f& a, const Vector3f& b ) { return _mm_cmpnle_ps(a.m128, b.m128); }
  __forceinline Vector3b ge_mask( const Vector3f& a, const Vector3f& b ) { return _mm_cmpnlt_ps(a.m128, b.m128); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Euclidian Space Operators
  ////////////////////////////////////////////////////////////////////////////////

#if defined(__SSE4_1__)
  __forceinline float dot ( const Vector3f& a, const Vector3f& b ) {
    return _mm_cvtss_f32(_mm_dp_ps(a,b,0x7F));
  }
#else
  __forceinline float dot ( const Vector3f& a, const Vector3f& b ) {
    return reduce_add(a*b);
  }
#endif

#if defined(__AVX2__)
  __forceinline Vector3f cross ( const Vector3f& a, const Vector3f& b ) 
  {
    ssef a0 = ssef(a);
    ssef b0 = shuffle<1,2,0,3>(ssef(b));
    ssef a1 = shuffle<1,2,0,3>(ssef(a));
    ssef b1 = ssef(b);
    return Vector3f(shuffle<1,2,0,3>(msub(a0,b0,a1*b1)));
  }
#else
  __forceinline Vector3f cross ( const Vector3f& a, const Vector3f& b ) 
  {
    ssef a0 = ssef(a);
    ssef b0 = shuffle<1,2,0,3>(ssef(b));
    ssef a1 = shuffle<1,2,0,3>(ssef(a));
    ssef b1 = ssef(b);
    return Vector3f(shuffle<1,2,0,3>(a0*b0-a1*b1));
  }
#endif

  __forceinline float   length   ( const Vector3f& a )                   { return sqrt(dot(a,a)); }
  __forceinline Vector3f normalize( const Vector3f& a )                   { return a*rsqrt(dot(a,a)); }
  __forceinline float   distance ( const Vector3f& a, const Vector3f& b ) { return length(a-b); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Select
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vector3f select( const Vector3b& s, const Vector3f& t, const Vector3f& f ) {
    return _mm_blendv_ps(f, t, s);
  }

  __forceinline const Vector3f select( bool s, const Vector3f& t, const Vector3f& f ) {
    __m128 mask = s ? _mm_castsi128_ps(_mm_cmpeq_epi32(_mm_setzero_si128(), _mm_setzero_si128())) : _mm_setzero_ps();
    return _mm_blendv_ps(f, t, mask);
  }

  __forceinline int maxDim ( const Vector3f& a ) 
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

  inline std::ostream& operator<<(std::ostream& cout, const Vector3f& a) {
    return cout << "(" << a.x << ", " << a.y << ", " << a.z << ")";
  }
}

#endif

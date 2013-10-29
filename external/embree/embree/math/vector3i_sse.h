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

#ifndef __EMBREE_VEC3I_SSE_H__
#define __EMBREE_VEC3I_SSE_H__

#include "simd/sse.h"
#include "math.h"

namespace embree
{
  ////////////////////////////////////////////////////////////////////////////////
  /// SSE Vector3i Type
  ////////////////////////////////////////////////////////////////////////////////

  struct __align(16) Vector3i
  {
    union {
      __m128i m128;
      struct { int x,y,z; int a; };
    };

    typedef int Scalar;
    enum { N = 3 };

    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vector3i( ) {}
    __forceinline Vector3i( const __m128i a ) : m128(a) {}
    __forceinline Vector3i( const Vector3i& other ) : m128(other.m128) {}
    __forceinline Vector3i& operator =(const Vector3i& other) { m128 = other.m128; return *this; }

    __forceinline Vector3i( const int a ) : m128(_mm_set1_epi32(a)) {}
    __forceinline explicit Vector3i( const int x, const int y, const int z) : m128(_mm_set_epi32(z, z, y, x)) {}
    __forceinline explicit Vector3i( const int* const a, const ssize_t stride = 1 ) : m128(_mm_set_epi32(0,a[2*stride],a[stride],a[0])) {}
    __forceinline explicit Vector3i( const __m128 a ) : m128(_mm_cvtps_epi32(a)) {}

    __forceinline operator const __m128i&( void ) const { return m128; }
    __forceinline operator       __m128i&( void )       { return m128; }

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vector3i( ZeroTy   ) : m128(_mm_setzero_si128()) {}
    __forceinline Vector3i( OneTy    ) : m128(_mm_set1_epi32(1)) {}
    __forceinline Vector3i( PosInfTy ) : m128(_mm_set1_epi32(pos_inf)) {}
    __forceinline Vector3i( NegInfTy ) : m128(_mm_set1_epi32(neg_inf)) {}

    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline const int& operator []( const size_t index ) const { assert(index < 3); return (&x)[index]; }
    __forceinline       int& operator []( const size_t index )       { assert(index < 3); return (&x)[index]; }
  };


  ////////////////////////////////////////////////////////////////////////////////
  /// Unary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vector3i operator +( const Vector3i& a ) { return a; }
  __forceinline const Vector3i operator -( const Vector3i& a ) { return _mm_sub_epi32(_mm_setzero_si128(), a.m128); }
  __forceinline const Vector3i abs       ( const Vector3i& a ) { return _mm_abs_epi32(a.m128); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vector3i operator +( const Vector3i& a, const Vector3i& b ) { return _mm_add_epi32(a.m128, b.m128); }
  __forceinline const Vector3i operator +( const Vector3i& a, const int        b ) { return a+Vector3i(b); }
  __forceinline const Vector3i operator +( const int        a, const Vector3i& b ) { return Vector3i(a)+b; }

  __forceinline const Vector3i operator -( const Vector3i& a, const Vector3i& b ) { return _mm_sub_epi32(a.m128, b.m128); }
  __forceinline const Vector3i operator -( const Vector3i& a, const int        b ) { return a-Vector3i(b); }
  __forceinline const Vector3i operator -( const int        a, const Vector3i& b ) { return Vector3i(a)-b; }

  __forceinline const Vector3i operator *( const Vector3i& a, const Vector3i& b ) { return _mm_mullo_epi32(a.m128, b.m128); }
  __forceinline const Vector3i operator *( const Vector3i& a, const int        b ) { return a * Vector3i(b); }
  __forceinline const Vector3i operator *( const int        a, const Vector3i& b ) { return Vector3i(a) * b; }

  __forceinline const Vector3i operator &( const Vector3i& a, const Vector3i& b ) { return _mm_and_si128(a.m128, b.m128); }
  __forceinline const Vector3i operator &( const Vector3i& a, const int        b ) { return a & Vector3i(b); }
  __forceinline const Vector3i operator &( const int        a, const Vector3i& b ) { return Vector3i(a) & b; }

  __forceinline const Vector3i operator |( const Vector3i& a, const Vector3i& b ) { return _mm_or_si128(a.m128, b.m128); }
  __forceinline const Vector3i operator |( const Vector3i& a, const int        b ) { return a | Vector3i(b); }
  __forceinline const Vector3i operator |( const int        a, const Vector3i& b ) { return Vector3i(a) | b; }

  __forceinline const Vector3i operator ^( const Vector3i& a, const Vector3i& b ) { return _mm_xor_si128(a.m128, b.m128); }
  __forceinline const Vector3i operator ^( const Vector3i& a, const int        b ) { return a ^ Vector3i(b); }
  __forceinline const Vector3i operator ^( const int        a, const Vector3i& b ) { return Vector3i(a) ^ b; }

  __forceinline const Vector3i operator <<( const Vector3i& a, const int n ) { return _mm_slli_epi32(a.m128, n); }
  __forceinline const Vector3i operator >>( const Vector3i& a, const int n ) { return _mm_srai_epi32(a.m128, n); }

  __forceinline const Vector3i sra ( const Vector3i& a, const int b ) { return _mm_srai_epi32(a.m128, b); }
  __forceinline const Vector3i srl ( const Vector3i& a, const int b ) { return _mm_srli_epi32(a.m128, b); }

  __forceinline const Vector3i min( const Vector3i& a, const Vector3i& b ) { return _mm_min_epi32(a.m128,b.m128); }
  __forceinline const Vector3i max( const Vector3i& a, const Vector3i& b ) { return _mm_max_epi32(a.m128,b.m128); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline Vector3i& operator +=( Vector3i& a, const Vector3i& b ) { return a = a + b; }
  __forceinline Vector3i& operator +=( Vector3i& a, const int32&     b ) { return a = a + b; }
  
  __forceinline Vector3i& operator -=( Vector3i& a, const Vector3i& b ) { return a = a - b; }
  __forceinline Vector3i& operator -=( Vector3i& a, const int32&     b ) { return a = a - b; }
  
  __forceinline Vector3i& operator *=( Vector3i& a, const Vector3i& b ) { return a = a * b; }
  __forceinline Vector3i& operator *=( Vector3i& a, const int32&     b ) { return a = a * b; }
  
  __forceinline Vector3i& operator &=( Vector3i& a, const Vector3i& b ) { return a = a & b; }
  __forceinline Vector3i& operator &=( Vector3i& a, const int32&     b ) { return a = a & b; }
  
  __forceinline Vector3i& operator |=( Vector3i& a, const Vector3i& b ) { return a = a | b; }
  __forceinline Vector3i& operator |=( Vector3i& a, const int32&     b ) { return a = a | b; }
  
  __forceinline Vector3i& operator <<=( Vector3i& a, const int32&    b ) { return a = a << b; }
  __forceinline Vector3i& operator >>=( Vector3i& a, const int32&    b ) { return a = a >> b; }

  ////////////////////////////////////////////////////////////////////////////////
  /// Reductions
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline int reduce_add(const Vector3i& v) { return v.x+v.y+v.z; }
  __forceinline int reduce_mul(const Vector3i& v) { return v.x*v.y*v.z; }
  __forceinline int reduce_min(const Vector3i& v) { return min(v.x,v.y,v.z); }
  __forceinline int reduce_max(const Vector3i& v) { return max(v.x,v.y,v.z); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline bool operator ==( const Vector3i& a, const Vector3i& b ) { return (_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpeq_epi32(a.m128, b.m128))) & 7) == 7; }
  __forceinline bool operator !=( const Vector3i& a, const Vector3i& b ) { return (_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpeq_epi32(a.m128, b.m128))) & 7) != 7; }
  __forceinline bool operator < ( const Vector3i& a, const Vector3i& b ) {
    if (a.x != b.x) return a.x < b.x;
    if (a.y != b.y) return a.y < b.y;
    if (a.z != b.z) return a.z < b.z;
    return false;
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Select
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vector3i select( const Vector3b& s, const Vector3i& t, const Vector3i& f ) {
    return _mm_castps_si128(_mm_blendv_ps(_mm_castsi128_ps(f), _mm_castsi128_ps(t), s));
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Output Operators
  ////////////////////////////////////////////////////////////////////////////////

  inline std::ostream& operator<<(std::ostream& cout, const Vector3i& a) {
    return cout << "(" << a.x << ", " << a.y << ", " << a.z << ")";
  }
}

#endif

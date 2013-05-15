// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#ifndef __EMBREE_VEC3B_SSE_H__
#define __EMBREE_VEC3B_SSE_H__

#include "simd/sse.h"
#include "math.h"

namespace embree
{
  ////////////////////////////////////////////////////////////////////////////////
  /// SSE Vec3b Type
  ////////////////////////////////////////////////////////////////////////////////

  template<> struct __align(16) Vec3<bool>
  {
    union {
      __m128 m128;
      struct { int x,y,z; int a; };
    };

    typedef int Scalar;
    enum { N = 3 };

    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vec3( ) {}
    __forceinline Vec3( const __m128  input ) : m128(input) {}
    __forceinline Vec3( const Vec3<bool>& other ) : m128(other.m128) {}
    __forceinline Vec3<bool>& operator =(const Vec3<bool>& other) { m128 = other.m128; return *this; }

    __forceinline explicit Vec3           ( bool  a )
      : m128(_mm_lookupmask_ps[(size_t(a) << 3) | (size_t(a) << 2) | (size_t(a) << 1) | size_t(a)]) {}
    __forceinline explicit Vec3           ( bool  a, bool  b, bool  c)
      : m128(_mm_lookupmask_ps[(size_t(c) << 2) | (size_t(b) << 1) | size_t(a)]) {}

    __forceinline operator const __m128&( void ) const { return m128; }
    __forceinline operator       __m128&( void )       { return m128; }

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vec3( FalseTy ) : m128(_mm_setzero_ps()) {}
    __forceinline Vec3( TrueTy  ) : m128(_mm_castsi128_ps(_mm_cmpeq_epi32(_mm_setzero_si128(), _mm_setzero_si128()))) {}

    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline const int& operator []( const size_t index ) const { assert(index < 3); return (&x)[index]; }
    __forceinline       int& operator []( const size_t index )       { assert(index < 3); return (&x)[index]; }
  };


  ////////////////////////////////////////////////////////////////////////////////
  /// Unary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vec3<bool> operator !( const Vec3<bool>& a ) { return _mm_xor_ps(a.m128, Vec3<bool>(True)); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vec3<bool> operator &( const Vec3<bool>& a, const Vec3<bool>& b ) { return _mm_and_ps(a.m128, b.m128); }
  __forceinline const Vec3<bool> operator |( const Vec3<bool>& a, const Vec3<bool>& b ) { return _mm_or_ps (a.m128, b.m128); }
  __forceinline const Vec3<bool> operator ^( const Vec3<bool>& a, const Vec3<bool>& b ) { return _mm_xor_ps(a.m128, b.m128); }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline const Vec3<bool> operator &=( Vec3<bool>& a, const Vec3<bool>& b ) { return a = a & b; }
  __forceinline const Vec3<bool> operator |=( Vec3<bool>& a, const Vec3<bool>& b ) { return a = a | b; }
  __forceinline const Vec3<bool> operator ^=( Vec3<bool>& a, const Vec3<bool>& b ) { return a = a ^ b; }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators + Select
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline bool operator ==( const Vec3<bool>& a, const Vec3<bool>& b ) { 
    return (_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpeq_epi32(_mm_castps_si128(a.m128), _mm_castps_si128(b.m128)))) & 7) == 7; 
  }
  __forceinline bool operator !=( const Vec3<bool>& a, const Vec3<bool>& b ) { 
    return (_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpeq_epi32(_mm_castps_si128(a.m128), _mm_castps_si128(b.m128)))) & 7) != 7; 
  }
  __forceinline bool operator < ( const Vec3<bool>& a, const Vec3<bool>& b ) {
    if (a.x != b.x) return a.x < b.x;
    if (a.y != b.y) return a.y < b.y;
    if (a.z != b.z) return a.z < b.z;
    return false;
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Output Operators
  ////////////////////////////////////////////////////////////////////////////////

  inline std::ostream& operator<<(std::ostream& cout, const Vec3<bool>& a) {
    return cout << "(" << (a.x ? "1" : "0") << ", " << (a.y ? "1" : "0") << ", " << (a.z ? "1" : "0") << ")";
  }
}

#endif

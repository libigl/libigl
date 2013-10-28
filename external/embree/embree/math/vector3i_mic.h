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

#ifndef __EMBREE_VEC3I_MIC_H__
#define __EMBREE_VEC3I_MIC_H__

#include "math.h"

namespace embree
{
  ////////////////////////////////////////////////////////////////////////////////
  /// MIC Vector3i Type
  ////////////////////////////////////////////////////////////////////////////////

  /*! 3-wide vectors emulated with 16-wide vectors. */
  struct __align(64) Vector3i 
  {
    typedef int Scalar;
    enum { N = 3 }; 
    union {
      __m512i m512; 
      struct { int x,y,z; int a; };
    };

    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vector3i( ) {}
    __forceinline Vector3i( const __m512i a ) : m512(a) {} 

    __forceinline Vector3i            ( const Vector3i& other ) { m512 = other.m512; }
    __forceinline Vector3i& operator =( const Vector3i& other ) { m512 = other.m512; return *this; }

    __forceinline Vector3i( const int a ) : m512(_mm512_set1_epi32(a)) {}
    __forceinline explicit Vector3i( const int x, const int y, const int z) : x(x), y(y), z(z) {}
    __forceinline explicit Vector3i( const int* const a, const ssize_t stride = 1 ) : x(a[0*stride]), y(a[1*stride]), z(a[2*stride]) {}

    __forceinline explicit Vector3i( const __m512 a ) : m512(_mm512_cvtfxpnt_round_adjustps_epi32(a,_MM_ROUND_MODE_NEAREST,_MM_EXPADJ_NONE)) {}

    __forceinline operator const __m512i&( void ) const { return m512; }
    __forceinline operator       __m512i&( void )       { return m512; }

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vector3i( ZeroTy   ) : m512(_mm512_setzero_epi32()) {}
    __forceinline Vector3i( OneTy    ) : m512(_mm512_set1_epi32(1)) {}
    __forceinline Vector3i( PosInfTy ) : m512(_mm512_set1_epi32(pos_inf)) {}
    __forceinline Vector3i( NegInfTy ) : m512(_mm512_set1_epi32(neg_inf)) {}

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
  __forceinline const Vector3i operator -( const Vector3i& a ) { return _mm512_sub_epi32(_mm512_setzero_epi32(), a.m512); }
  //__forceinline const Vector3i abs       ( const Vector3i& a ) { return _mm512_abs_epi32(a.m512); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vector3i operator +( const Vector3i& a, const Vector3i& b ) { return _mm512_add_epi32(a.m512, b.m512); }
  __forceinline const Vector3i operator +( const Vector3i& a, const int      b ) { return a+Vector3i(b); }
  __forceinline const Vector3i operator +( const int      a, const Vector3i& b ) { return Vector3i(a)+b; }

  __forceinline const Vector3i operator -( const Vector3i& a, const Vector3i& b ) { return _mm512_sub_epi32(a.m512, b.m512); }
  __forceinline const Vector3i operator -( const Vector3i& a, const int      b ) { return a-Vector3i(b); }
  __forceinline const Vector3i operator -( const int      a, const Vector3i& b ) { return Vector3i(a)-b; }

  __forceinline const Vector3i operator *( const Vector3i& a, const Vector3i& b ) { return _mm512_mullo_epi32(a.m512, b.m512); }
  __forceinline const Vector3i operator *( const Vector3i& a, const int      b ) { return a * Vector3i(b); }
  __forceinline const Vector3i operator *( const int      a, const Vector3i& b ) { return Vector3i(a) * b; }

  __forceinline const Vector3i operator &( const Vector3i& a, const Vector3i& b ) { return _mm512_and_epi32(a.m512, b.m512); }
  __forceinline const Vector3i operator &( const Vector3i& a, const int      b ) { return a & Vector3i(b); }
  __forceinline const Vector3i operator &( const int      a, const Vector3i& b ) { return Vector3i(a) & b; }

  __forceinline const Vector3i operator |( const Vector3i& a, const Vector3i& b ) { return _mm512_or_epi32(a.m512, b.m512); }
  __forceinline const Vector3i operator |( const Vector3i& a, const int      b ) { return a | Vector3i(b); }
  __forceinline const Vector3i operator |( const int      a, const Vector3i& b ) { return Vector3i(a) | b; }

  __forceinline const Vector3i operator ^( const Vector3i& a, const Vector3i& b ) { return _mm512_xor_epi32(a.m512, b.m512); }
  __forceinline const Vector3i operator ^( const Vector3i& a, const int      b ) { return a ^ Vector3i(b); }
  __forceinline const Vector3i operator ^( const int      a, const Vector3i& b ) { return Vector3i(a) ^ b; }

  __forceinline const Vector3i operator <<( const Vector3i& a, const int n ) { return _mm512_slli_epi32(a.m512, n); }
  __forceinline const Vector3i operator >>( const Vector3i& a, const int n ) { return _mm512_srai_epi32(a.m512, n); }

  __forceinline const Vector3i sra ( const Vector3i& a, const int b ) { return _mm512_srai_epi32(a.m512, b); }
  __forceinline const Vector3i srl ( const Vector3i& a, const int b ) { return _mm512_srli_epi32(a.m512, b); }

  __forceinline const Vector3i min( const Vector3i& a, const Vector3i& b ) { return _mm512_min_epi32(a.m512,b.m512); }
  __forceinline const Vector3i max( const Vector3i& a, const Vector3i& b ) { return _mm512_max_epi32(a.m512,b.m512); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline Vector3i& operator +=( Vector3i& a, const Vector3i& b ) { return a = a + b; }
  __forceinline Vector3i& operator +=( Vector3i& a, const int32&   b ) { return a = a + b; }
  
  __forceinline Vector3i& operator -=( Vector3i& a, const Vector3i& b ) { return a = a - b; }
  __forceinline Vector3i& operator -=( Vector3i& a, const int32&   b ) { return a = a - b; }
  
  __forceinline Vector3i& operator *=( Vector3i& a, const Vector3i& b ) { return a = a * b; }
  __forceinline Vector3i& operator *=( Vector3i& a, const int32&   b ) { return a = a * b; }
  
  __forceinline Vector3i& operator &=( Vector3i& a, const Vector3i& b ) { return a = a & b; }
  __forceinline Vector3i& operator &=( Vector3i& a, const int32&   b ) { return a = a & b; }
  
  __forceinline Vector3i& operator |=( Vector3i& a, const Vector3i& b ) { return a = a | b; }
  __forceinline Vector3i& operator |=( Vector3i& a, const int32&   b ) { return a = a | b; }
  
  __forceinline Vector3i& operator <<=( Vector3i& a, const int32&  b ) { return a = a << b; }
  __forceinline Vector3i& operator >>=( Vector3i& a, const int32&  b ) { return a = a >> b; }

  ////////////////////////////////////////////////////////////////////////////////
  /// Reductions
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline int reduce_add(const Vector3i& a) { 
    __m512i b = _mm512_mask_blend_epi32(0x8,a.m512,_mm512_set1_epi32(0));
    __m512i c = _mm512_add_epi32(b, _mm512_swizzle_epi32(b, _MM_SWIZ_REG_BADC));
    __m512i d = _mm512_add_epi32(c, _mm512_swizzle_epi32(c, _MM_SWIZ_REG_CDAB));
    return *(int*)&d;
  } 

  __forceinline int reduce_mul(const Vector3i& a) { 
    __m512i b = _mm512_mask_blend_epi32(0x8,a.m512,_mm512_set1_epi32(1));
    __m512i c = _mm512_mullo_epi32(b, _mm512_swizzle_epi32(b, _MM_SWIZ_REG_BADC));
    __m512i d = _mm512_mullo_epi32(c, _mm512_swizzle_epi32(c, _MM_SWIZ_REG_CDAB));
    return *(int*)&d;
  }

  __forceinline int reduce_min(const Vector3i& a) { 
    __m512i b = _mm512_mask_blend_epi32(0x8,a.m512,_mm512_set1_epi32(pos_inf));
    __m512i c = _mm512_min_epi32(b, _mm512_swizzle_epi32(b, _MM_SWIZ_REG_BADC));
    __m512i d = _mm512_min_epi32(c, _mm512_swizzle_epi32(c, _MM_SWIZ_REG_CDAB));
    return *(int*)&d;
  } 

  __forceinline int reduce_max(const Vector3i& a) { 
    __m512i b = _mm512_mask_blend_epi32(0x8,a.m512,_mm512_set1_epi32(neg_inf));
    __m512i c = _mm512_max_epi32(b, _mm512_swizzle_epi32(b, _MM_SWIZ_REG_BADC));
    __m512i d = _mm512_max_epi32(c, _mm512_swizzle_epi32(c, _MM_SWIZ_REG_CDAB));
    return *(int*)&d;
  } 

  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline bool operator ==( const Vector3i& a, const Vector3i& b ) { return (_mm512_cmp_epi32_mask(a,b,_MM_CMPINT_EQ) & 7) == 7; }
  __forceinline bool operator !=( const Vector3i& a, const Vector3i& b ) { return (_mm512_cmp_epi32_mask(a,b,_MM_CMPINT_NE) & 7) != 0; }
  __forceinline bool operator < ( const Vector3i& a, const Vector3i& b ) {
    if (a.x != b.x) return a.x < b.x;
    if (a.y != b.y) return a.y < b.y;
    if (a.z != b.z) return a.z < b.z;
    return false;
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Select
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vector3i select( bool s, const Vector3i& t, const Vector3i& f ) {
    return _mm512_mask_blend_epi32(s ? 0xF : 0x0, f, t);
  }

  __forceinline const Vector3i select( const Vector3b& s, const Vector3i& t, const Vector3i& f ) {
    return _mm512_mask_blend_epi32(s, f, t);
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Output Operators
  ////////////////////////////////////////////////////////////////////////////////

  inline std::ostream& operator<<(std::ostream& cout, const Vector3i& a) {
    return cout << "(" << a.x << ", " << a.y << ", " << a.z << ")";
  }
}

#endif

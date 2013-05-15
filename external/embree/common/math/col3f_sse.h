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

#ifndef __EMBREE_COL3F_SSE_H__
#define __EMBREE_COL3F_SSE_H__

#include "../simd/sse.h"
#include "../simd/sse_special.h"

namespace embree
{
   ////////////////////////////////////////////////////////////////////////////////
  /// SSE RGB Color Class
  ////////////////////////////////////////////////////////////////////////////////

  template<> struct Col3<float>
  {
    union {
      __m128 m128;
      struct { float r,g,b; int dummy; };
    };

    ////////////////////////////////////////////////////////////////////////////////
    /// Construction
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Col3           ( )                   { }
    __forceinline Col3           ( const __m128 a           ) : m128(a) {}
    __forceinline Col3           ( const Col3<float>& other ) : m128(other.m128) {}
    __forceinline Col3<float>& operator=( const Col3<float>& other ) { m128 = other.m128; return *this; }

    __forceinline explicit Col3 (const float v)                               : m128(_mm_set1_ps(v)) {}
    __forceinline explicit Col3 (const float* v, int stride = 1)              : m128(_mm_set_ps(0.0f,v[2*stride],v[stride],v[0])) {}
    __forceinline          Col3 (const float r, const float g, const float b) : m128(_mm_set_ps(0.0f,b,g,r)) {}

    __forceinline operator const __m128&( void ) const { return m128; }
    __forceinline operator       __m128&( void )       { return m128; }

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Col3( ZeroTy   ) : m128(_mm_set1_ps(0.0f)) {}
    __forceinline Col3( OneTy    ) : m128(_mm_set1_ps(1.0f)) {}
    __forceinline Col3( PosInfTy ) : m128(_mm_set1_ps(pos_inf)) {}
    __forceinline Col3( NegInfTy ) : m128(_mm_set1_ps(neg_inf)) {}
  };

  ////////////////////////////////////////////////////////////////////////////////
  /// Unary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Col3<float> operator +( const Col3<float>& a ) { return a; }
  __forceinline const Col3<float> operator -( const Col3<float>& a ) {
    const __m128 mask = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
    return _mm_xor_ps(a.m128, mask);
  }
  __forceinline const Col3<float> abs  ( const Col3<float>& a ) {
    const __m128 mask = _mm_castsi128_ps(_mm_set1_epi32(0x7fffffff));
    return _mm_and_ps(a.m128, mask);
  }
  __forceinline const Col3<float> rcp  ( const Col3<float>& a ) {
    const Col3<float> r = _mm_rcp_ps(a.m128);
    return _mm_sub_ps(_mm_add_ps(r, r), _mm_mul_ps(_mm_mul_ps(r, r), a));
  }
  __forceinline const Col3<float> rsqrt( const Col3<float>& a ) {
    __m128 r = _mm_rsqrt_ps(a.m128);
    return _mm_add_ps(_mm_mul_ps(_mm_set1_ps(1.5f),r), _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(a, _mm_set1_ps(-0.5f)), r), _mm_mul_ps(r, r)));
  }
  __forceinline const Col3<float> sqrt ( const Col3<float>& a ) { return _mm_sqrt_ps(a.m128); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Col3<float> operator +( const Col3<float>& a, const Col3<float>& b ) { return _mm_add_ps(a.m128, b.m128); }
  __forceinline const Col3<float> operator -( const Col3<float>& a, const Col3<float>& b ) { return _mm_sub_ps(a.m128, b.m128); }
  __forceinline const Col3<float> operator *( const Col3<float>& a, const Col3<float>& b ) { return _mm_mul_ps(a.m128, b.m128); }
  __forceinline const Col3<float> operator *( const Col3<float>& a, const float b        ) { return a * Col3<float>(b); }
  __forceinline const Col3<float> operator *( const float a       , const Col3<float>& b ) { return Col3<float>(a) * b; }
  __forceinline const Col3<float> operator /( const Col3<float>& a, const Col3<float>& b ) { return a * rcp(b); }
  __forceinline const Col3<float> operator /( const Col3<float>& a, const float b ) { return a * rcp(b); }

  __forceinline const Col3<float> min( const Col3<float>& a, const Col3<float>& b ) { return _mm_min_ps(a.m128,b.m128); }
  __forceinline const Col3<float> max( const Col3<float>& a, const Col3<float>& b ) { return _mm_max_ps(a.m128,b.m128); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Col3<float> operator+=(Col3<float>& a, const Col3<float>& b) { return a = a + b; }
  __forceinline const Col3<float> operator-=(Col3<float>& a, const Col3<float>& b) { return a = a - b; }
  __forceinline const Col3<float> operator*=(Col3<float>& a, const Col3<float>& b) { return a = a * b; }
  __forceinline const Col3<float> operator/=(Col3<float>& a, const Col3<float>& b) { return a = a / b; }
  __forceinline const Col3<float> operator*=(Col3<float>& a, const float b      ) { return a = a * b; }
  __forceinline const Col3<float> operator/=(Col3<float>& a, const float b      ) { return a = a / b; }

  ////////////////////////////////////////////////////////////////////////////////
  /// Reductions
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline float reduce_add(const Col3<float>& v) { return v.r+v.g+v.b; }
  __forceinline float reduce_mul(const Col3<float>& v) { return v.r*v.g*v.b; }
  __forceinline float reduce_min(const Col3<float>& v) { return min(v.r,v.g,v.b); }
  __forceinline float reduce_max(const Col3<float>& v) { return max(v.r,v.g,v.b); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline bool operator ==( const Col3<float>& a, const Col3<float>& b ) { return (_mm_movemask_ps(_mm_cmpeq_ps (a.m128, b.m128)) & 7) == 7; }
  __forceinline bool operator !=( const Col3<float>& a, const Col3<float>& b ) { return (_mm_movemask_ps(_mm_cmpneq_ps(a.m128, b.m128)) & 7) != 0; }
  __forceinline bool operator < ( const Col3<float>& a, const Col3<float>& b ) {
    if (a.r != b.r) return a.r < b.r;
    if (a.g != b.g) return a.g < b.g;
    if (a.b != b.b) return a.b < b.b;
    return false;
  }

   ////////////////////////////////////////////////////////////////////////////////
  /// Select
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Col3<float> select( bool s, const Col3<float>& t, const Col3<float>& f ) {
    __m128 mask = s ? _mm_castsi128_ps(_mm_cmpeq_epi32(_mm_setzero_si128(), _mm_setzero_si128())) : _mm_setzero_ps();
    return _mm_blendv_ps(f, t, mask);
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Special Operators
  ////////////////////////////////////////////////////////////////////////////////

  /*! computes luminance of a color */
  __forceinline float luminance (const Col3<float>& a) { return 0.212671f*a.r + 0.715160f*a.g + 0.072169f*a.b; }

  __forceinline Col3<float> exp (const Col3<float>& a) { return Col3<float>(exp_ps(a.m128)); }
  __forceinline Col3<float> log (const Col3<float>& a) { return Col3<float>(log_ps(a.m128)); }

  /*! output operator */
  inline std::ostream& operator<<(std::ostream& cout, const Col3<float>& a) {
    return cout << "(" << a.r << ", " << a.g << ", " << a.b << ")";
  }
}

#endif

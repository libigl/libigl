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

#ifndef __EMBREE_VEC3F_SSE_H__
#define __EMBREE_VEC3F_SSE_H__

#include "simd/sse.h"
#include "math.h"

namespace embree
{
  ////////////////////////////////////////////////////////////////////////////////
  /// SSE Vec3f Type
  ////////////////////////////////////////////////////////////////////////////////

  template<> struct __align(16) Vec3<float>
  {
    union {
      __m128 m128;
      struct { float x,y,z; int a; };
    };

    typedef float Scalar;
    enum { N = 3 };

    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vec3( ) {}
    __forceinline Vec3( const __m128 a ) : m128(a) {}
    __forceinline Vec3( const Vec3<float>& other ) : m128(other.m128) {}
    __forceinline Vec3<float>& operator =(const Vec3<float>& other) { m128 = other.m128; return *this; }

    __forceinline explicit Vec3( const float a ) : m128(_mm_set1_ps(a)) {}
    __forceinline explicit Vec3( const float x, const float y, const float z) : m128(_mm_set_ps(z, z, y, x)) {}
    __forceinline explicit Vec3( const float* const a, const ssize_t stride = 1 ) : m128(_mm_set_ps(0.0f,a[2*stride],a[stride],a[0])) {}
    __forceinline explicit Vec3( const __m128i a ) : m128(_mm_cvtepi32_ps(a)) {}

    __forceinline operator const __m128&( void ) const { return m128; }
    __forceinline operator       __m128&( void )       { return m128; }

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vec3( ZeroTy   ) : m128(_mm_setzero_ps()) {}
    __forceinline Vec3( OneTy    ) : m128(_mm_set1_ps(1.0f)) {}
    __forceinline Vec3( PosInfTy ) : m128(_mm_set1_ps(pos_inf)) {}
    __forceinline Vec3( NegInfTy ) : m128(_mm_set1_ps(neg_inf)) {}

    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline const float& operator []( const size_t index ) const { assert(index < 3); return (&x)[index]; }
    __forceinline       float& operator []( const size_t index )       { assert(index < 3); return (&x)[index]; }
  };


  ////////////////////////////////////////////////////////////////////////////////
  /// Unary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vec3<float> operator +( const Vec3<float>& a ) { return a; }
  __forceinline const Vec3<float> operator -( const Vec3<float>& a ) {
    const __m128 mask = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
    return _mm_xor_ps(a.m128, mask);
  }
  __forceinline const Vec3<float> abs  ( const Vec3<float>& a ) {
    const __m128 mask = _mm_castsi128_ps(_mm_set1_epi32(0x7fffffff));
    return _mm_and_ps(a.m128, mask);
  }
  __forceinline const Vec3<float> sign ( const Vec3<float>& a ) {
    return _mm_blendv_ps(Vec3<float>(one), -Vec3<float>(one), _mm_cmplt_ps (a,Vec3<float>(zero))); 
  }
  __forceinline const Vec3<float> rcp  ( const Vec3<float>& a ) {
    const Vec3<float> r = _mm_rcp_ps(a.m128);
    return _mm_sub_ps(_mm_add_ps(r, r), _mm_mul_ps(_mm_mul_ps(r, r), a));
  }
  __forceinline const Vec3<float> sqrt ( const Vec3<float>& a ) { return _mm_sqrt_ps(a.m128); }
  __forceinline const Vec3<float> sqr  ( const Vec3<float>& a ) { return _mm_mul_ps(a,a); }
  __forceinline const Vec3<float> rsqrt( const Vec3<float>& a ) {
    __m128 r = _mm_rsqrt_ps(a.m128);
    return _mm_add_ps(_mm_mul_ps(_mm_set1_ps(1.5f),r), _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(a, _mm_set1_ps(-0.5f)), r), _mm_mul_ps(r, r)));
  }
  __forceinline const Vec3<float> rcp_safe( const Vec3<float>& a ) 
  {
    __m128 mask = _mm_cmpeq_ps (a.m128, _mm_setzero_ps());

	/*! workaround for compiler bug in VS2008 */
#if defined(_MSC_VER) && (_MSC_VER < 1600)
    return _mm_or_ps(_mm_and_ps(mask, _mm_set1_ps(std::numeric_limits<float>::max())), _mm_andnot_ps(mask, rcp(a))); 
#else
	return _mm_blendv_ps(rcp(a), _mm_set1_ps(std::numeric_limits<float>::max()), mask);
#endif
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vec3<float> operator +( const Vec3<float>& a, const Vec3<float>& b ) { return _mm_add_ps(a.m128, b.m128); }
  __forceinline const Vec3<float> operator -( const Vec3<float>& a, const Vec3<float>& b ) { return _mm_sub_ps(a.m128, b.m128); }
  __forceinline const Vec3<float> operator *( const Vec3<float>& a, const Vec3<float>& b ) { return _mm_mul_ps(a.m128, b.m128); }
  __forceinline const Vec3<float> operator *( const Vec3<float>& a, const float b ) { return a * Vec3<float>(b); }
  __forceinline const Vec3<float> operator *( const float a, const Vec3<float>& b ) { return Vec3<float>(a) * b; }
  __forceinline const Vec3<float> operator /( const Vec3<float>& a, const Vec3<float>& b ) { return _mm_div_ps(a.m128,b.m128); }
  __forceinline const Vec3<float> operator /( const Vec3<float>& a, const float b        ) { return _mm_div_ps(a.m128,_mm_set1_ps(b)); }
  __forceinline const Vec3<float> operator /( const        float a, const Vec3<float>& b ) { return _mm_div_ps(_mm_set1_ps(a),b.m128); }

  __forceinline const Vec3<float> min( const Vec3<float>& a, const Vec3<float>& b ) { return _mm_min_ps(a.m128,b.m128); }
  __forceinline const Vec3<float> max( const Vec3<float>& a, const Vec3<float>& b ) { return _mm_max_ps(a.m128,b.m128); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline Vec3<float>& operator +=( Vec3<float>& a, const Vec3<float>& b ) { return a = a + b; }
  __forceinline Vec3<float>& operator -=( Vec3<float>& a, const Vec3<float>& b ) { return a = a - b; }
  __forceinline Vec3<float>& operator *=( Vec3<float>& a, const float b ) { return a = a * b; }
  __forceinline Vec3<float>& operator /=( Vec3<float>& a, const float b ) { return a = a / b; }

  ////////////////////////////////////////////////////////////////////////////////
  /// Reductions
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline float reduce_add(const Vec3<float>& v) { return v.x+v.y+v.z; }
  __forceinline float reduce_mul(const Vec3<float>& v) { return v.x*v.y*v.z; }
  __forceinline float reduce_min(const Vec3<float>& v) { return min(v.x,v.y,v.z); }
  __forceinline float reduce_max(const Vec3<float>& v) { return max(v.x,v.y,v.z); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline bool operator ==( const Vec3<float>& a, const Vec3<float>& b ) { return (_mm_movemask_ps(_mm_cmpeq_ps (a.m128, b.m128)) & 7) == 7; }
  __forceinline bool operator !=( const Vec3<float>& a, const Vec3<float>& b ) { return (_mm_movemask_ps(_mm_cmpneq_ps(a.m128, b.m128)) & 7) != 0; }
  __forceinline bool operator < ( const Vec3<float>& a, const Vec3<float>& b ) {
    if (a.x != b.x) return a.x < b.x;
    if (a.y != b.y) return a.y < b.y;
    if (a.z != b.z) return a.z < b.z;
    return false;
  }

  __forceinline Vec3<bool> eq_mask( const Vec3<float>& a, const Vec3<float>& b ) { return _mm_cmpeq_ps (a.m128, b.m128); }
  __forceinline Vec3<bool> neq_mask(const Vec3<float>& a, const Vec3<float>& b ) { return _mm_cmpneq_ps(a.m128, b.m128); }
  __forceinline Vec3<bool> lt_mask( const Vec3<float>& a, const Vec3<float>& b ) { return _mm_cmplt_ps (a.m128, b.m128); }
  __forceinline Vec3<bool> le_mask( const Vec3<float>& a, const Vec3<float>& b ) { return _mm_cmple_ps (a.m128, b.m128); }
  __forceinline Vec3<bool> gt_mask( const Vec3<float>& a, const Vec3<float>& b ) { return _mm_cmpnle_ps(a.m128, b.m128); }
  __forceinline Vec3<bool> ge_mask( const Vec3<float>& a, const Vec3<float>& b ) { return _mm_cmpnlt_ps(a.m128, b.m128); }

 ////////////////////////////////////////////////////////////////////////////////
  /// Euclidian Space Operators
  ////////////////////////////////////////////////////////////////////////////////

#if defined(__SSE4_1__)
  __forceinline float dot ( const Vec3<float>& a, const Vec3<float>& b ) {
    return _mm_cvtss_f32(_mm_dp_ps(a,b,0x7F));
  }
#else
  __forceinline float dot ( const Vec3<float>& a, const Vec3<float>& b ) {
    return reduce_add(a*b);
  }
#endif

  __forceinline Vec3<float> cross ( const Vec3<float>& a, const Vec3<float>& b ) {
    ssef a0 = shuffle<1,2,0,3>(ssef(a));
    ssef b0 = shuffle<2,0,1,3>(ssef(b));
    ssef a1 = shuffle<2,0,1,3>(ssef(a));
    ssef b1 = shuffle<1,2,0,3>(ssef(b));
    return Vec3<float>(a0*b0-a1*b1);
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Select
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vec3<float> select( bool s, const Vec3<float>& t, const Vec3<float>& f ) {
    __m128 mask = s ? _mm_castsi128_ps(_mm_cmpeq_epi32(_mm_setzero_si128(), _mm_setzero_si128())) : _mm_setzero_ps();
    return _mm_blendv_ps(f, t, mask);
  }

  __forceinline int maxDim ( const Vec3<float>& a ) 
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

  inline std::ostream& operator<<(std::ostream& cout, const Vec3<float>& a) {
    return cout << "(" << a.x << ", " << a.y << ", " << a.z << ")";
  }
}

#endif

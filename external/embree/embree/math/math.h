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

#ifndef __EMBREE_MATH_H__
#define __EMBREE_MATH_H__

#include "sys/platform.h"
#include "sys/intrinsics.h"

#include <cmath>
#include <float.h>
#include <emmintrin.h>
#include <xmmintrin.h>

namespace embree
{
#if defined(__WIN32__)
  __forceinline bool finite ( const float x ) { return _finite(x) != 0; }
#endif

  __forceinline float sign ( const float x ) { return x<0?-1.0f:1.0f; }
  __forceinline float sqr  ( const float x ) { return x*x; }

#ifndef __MIC__
  __forceinline float rcp  ( const float x ) 
  { 
    const __m128 vx = _mm_set_ss(x);
    const __m128 r = _mm_rcp_ps(vx);
    return _mm_cvtss_f32(_mm_sub_ps(_mm_add_ps(r, r), _mm_mul_ps(_mm_mul_ps(r, r), vx)));
  }
  __forceinline float signmsk ( const float x ) { 
    return _mm_cvtss_f32(_mm_and_ps(_mm_set_ss(x),_mm_castsi128_ps(_mm_set1_epi32(0x80000000))));
  }
  __forceinline float xorf( const float x, const float y ) { 
    return _mm_cvtss_f32(_mm_xor_ps(_mm_set_ss(x),_mm_set_ss(y)));
  }
  __forceinline float andf( const float x, const unsigned y ) { 
    return _mm_cvtss_f32(_mm_and_ps(_mm_set_ss(x),_mm_castsi128_ps(_mm_set1_epi32(y))));
  }
  __forceinline float rsqrt( const float x ) { 
    const __m128 a = _mm_set_ss(x);
    const __m128 r = _mm_rsqrt_ps(a);
    const __m128 c = _mm_add_ps(_mm_mul_ps(_mm_set_ps(1.5f, 1.5f, 1.5f, 1.5f), r),
                                _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(a, _mm_set_ps(-0.5f, -0.5f, -0.5f, -0.5f)), r), _mm_mul_ps(r, r)));
    return _mm_cvtss_f32(c);
  }
#else
  __forceinline float rcp  ( const float x ) { return 1.0f/x; }
  __forceinline float signmsk ( const float x ) { return cast_i2f(cast_f2i(x)&0x80000000); }
  __forceinline float xorf( const float x, const float y ) { return cast_i2f(cast_f2i(x) ^ cast_f2i(y)); }
  __forceinline float andf( const float x, const float y ) { return cast_i2f(cast_f2i(x) & cast_f2i(y)); }
  __forceinline float rsqrt( const float x ) { return 1.0f/sqrtf(x); }
#endif

  __forceinline float abs  ( const float x ) { return ::fabsf(x); }
  __forceinline float acos ( const float x ) { return ::acosf (x); }
  __forceinline float asin ( const float x ) { return ::asinf (x); }
  __forceinline float atan ( const float x ) { return ::atanf (x); }
  __forceinline float atan2( const float y, const float x ) { return ::atan2f(y, x); }
  __forceinline float cos  ( const float x ) { return ::cosf  (x); }
  __forceinline float cosh ( const float x ) { return ::coshf (x); }
  __forceinline float exp  ( const float x ) { return ::expf  (x); }
  __forceinline float fmod ( const float x, const float y ) { return ::fmodf (x, y); }
  __forceinline float log  ( const float x ) { return ::logf  (x); }
  __forceinline float log10( const float x ) { return ::log10f(x); }
  __forceinline float pow  ( const float x, const float y ) { return ::powf  (x, y); }
  __forceinline float sin  ( const float x ) { return ::sinf  (x); }
  __forceinline float sinh ( const float x ) { return ::sinhf (x); }
  __forceinline float sqrt ( const float x ) { return ::sqrtf (x); }
  __forceinline float tan  ( const float x ) { return ::tanf  (x); }
  __forceinline float tanh ( const float x ) { return ::tanhf (x); }
  __forceinline float floor( const float x ) { return ::floorf (x); }
  __forceinline float ceil ( const float x ) { return ::ceilf (x); }

  __forceinline double abs  ( const double x ) { return ::fabs(x); }
  __forceinline double sign ( const double x ) { return x<0?-1.0:1.0; }
  __forceinline double acos ( const double x ) { return ::acos (x); }
  __forceinline double asin ( const double x ) { return ::asin (x); }
  __forceinline double atan ( const double x ) { return ::atan (x); }
  __forceinline double atan2( const double y, const double x ) { return ::atan2(y, x); }
  __forceinline double cos  ( const double x ) { return ::cos  (x); }
  __forceinline double cosh ( const double x ) { return ::cosh (x); }
  __forceinline double exp  ( const double x ) { return ::exp  (x); }
  __forceinline double fmod ( const double x, const double y ) { return ::fmod (x, y); }
  __forceinline double log  ( const double x ) { return ::log  (x); }
  __forceinline double log10( const double x ) { return ::log10(x); }
  __forceinline double pow  ( const double x, const double y ) { return ::pow  (x, y); }
  __forceinline double rcp  ( const double x ) { return 1.0/x; }
  __forceinline double rsqrt( const double x ) { return 1.0/::sqrt(x); }
  __forceinline double sin  ( const double x ) { return ::sin  (x); }
  __forceinline double sinh ( const double x ) { return ::sinh (x); }
  __forceinline double sqr  ( const double x ) { return x*x; }
  __forceinline double sqrt ( const double x ) { return ::sqrt (x); }
  __forceinline double tan  ( const double x ) { return ::tan  (x); }
  __forceinline double tanh ( const double x ) { return ::tanh (x); }
  __forceinline double floor( const double x ) { return ::floor (x); }
  __forceinline double ceil ( const double x ) { return ::ceil (x); }

  template<typename T> __forceinline T min(const T& a, const T& b)                                     { return a<b? a:b; }
  template<typename T> __forceinline T min(const T& a, const T& b, const T& c)                         { return min(min(a,b),c); }
  template<typename T> __forceinline T min(const T& a, const T& b, const T& c, const T& d)             { return min(min(a,b),min(c,d)); }
  template<typename T> __forceinline T min(const T& a, const T& b, const T& c, const T& d, const T& e) { return min(min(min(a,b),min(c,d)),e); }

  template<typename T> __forceinline T max(const T& a, const T& b)                                     { return a<b? b:a; }
  template<typename T> __forceinline T max(const T& a, const T& b, const T& c)                         { return max(max(a,b),c); }
  template<typename T> __forceinline T max(const T& a, const T& b, const T& c, const T& d)             { return max(max(a,b),max(c,d)); }
  template<typename T> __forceinline T max(const T& a, const T& b, const T& c, const T& d, const T& e) { return max(max(max(a,b),max(c,d)),e); }

  template<typename T> __forceinline T clamp(const T& x, const T& lower = T(zero), const T& upper = T(one)) { return max(lower, min(x,upper)); }

  template<typename T> __forceinline T  deg2rad ( const T& x )  { return x * T(1.74532925199432957692e-2f); }
  template<typename T> __forceinline T  rad2deg ( const T& x )  { return x * T(5.72957795130823208768e1f); }
  template<typename T> __forceinline T  sin2cos ( const T& x )  { return sqrt(max(T(zero),T(one)-x*x)); }
  template<typename T> __forceinline T  cos2sin ( const T& x )  { return sin2cos(x); }

#if defined(__AVX2__)
  __forceinline float madd  ( const float a, const float b, const float c) { return _mm_cvtss_f32(_mm_fmadd_ss(_mm_set_ss(a),_mm_set_ss(b),_mm_set_ss(c))); }
  __forceinline float msub  ( const float a, const float b, const float c) { return _mm_cvtss_f32(_mm_fmsub_ss(_mm_set_ss(a),_mm_set_ss(b),_mm_set_ss(c))); }
  __forceinline float nmadd ( const float a, const float b, const float c) { return _mm_cvtss_f32(_mm_fnmadd_ss(_mm_set_ss(a),_mm_set_ss(b),_mm_set_ss(c))); }
  __forceinline float nmsub ( const float a, const float b, const float c) { return _mm_cvtss_f32(_mm_fnmsub_ss(_mm_set_ss(a),_mm_set_ss(b),_mm_set_ss(c))); }
#else
  __forceinline float madd  ( const float a, const float b, const float c) { return a*b+c; }
  __forceinline float msub  ( const float a, const float b, const float c) { return a*b-c; }
  __forceinline float nmadd ( const float a, const float b, const float c) { return -a*b-c;}
  __forceinline float nmsub ( const float a, const float b, const float c) { return c-a*b; }
#endif

  /*! random functions */
  template<typename T> T          random() { return T(0); }
  template<> __forceinline int    random() { return int(rand()); }
  template<> __forceinline uint32 random() { return uint32(rand()); }
  template<> __forceinline float  random() { return random<uint32>()/float(RAND_MAX); }
  template<> __forceinline double random() { return random<uint32>()/double(RAND_MAX); }

  /*! selects */
  __forceinline bool  select(bool s, bool  t , bool f) { return s ? t : f; }
  __forceinline int   select(bool s, int   t,   int f) { return s ? t : f; }
  __forceinline float select(bool s, float t, float f) { return s ? t : f; }
}

#endif

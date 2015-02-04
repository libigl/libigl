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

#include "sys/platform.h"
#include "math/math.h"

#if defined __SSE__
#include "simd/sse.h"
#endif

#if defined __AVX__
#include "simd/avx.h"
#endif

//#if defined __MIC__
//#include "simd/sse_mic.h"
//#endif

#if defined __MIC__
#include "simd/mic.h"
#endif

namespace embree
{
  struct Vec3fa;

  ////////////////////////////////////////////////////////////////////////////////
  /// Generic 3D vector Class
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> struct Vec3
  {
    T x, y, z;

    typedef T Scalar;
    enum { N  = 3 };

    ////////////////////////////////////////////////////////////////////////////////
    /// Construction
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vec3 ( ) {}
    __forceinline Vec3 ( const T& a                         ) : x(a), y(a), z(a) {}
    __forceinline Vec3 ( const T& x, const T& y, const T& z ) : x(x), y(y), z(z) {}

    __forceinline Vec3     ( const Vec3& other ) { x = other.x; y = other.y; z = other.z; }
    __forceinline Vec3     ( const Vec3fa& other );

    template<typename T1> __forceinline Vec3( const Vec3<T1>& a ) : x(T(a.x)), y(T(a.y)), z(T(a.z)) {}
    template<typename T1> __forceinline Vec3& operator =(const Vec3<T1>& other) { x = other.x; y = other.y; z = other.z; return *this; }

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vec3( ZeroTy   ) : x(zero), y(zero), z(zero) {}
    __forceinline Vec3( OneTy    ) : x(one),  y(one),  z(one) {}
    __forceinline Vec3( PosInfTy ) : x(pos_inf), y(pos_inf), z(pos_inf) {}
    __forceinline Vec3( NegInfTy ) : x(neg_inf), y(neg_inf), z(neg_inf) {}

    __forceinline const T& operator []( const size_t axis ) const { assert(axis < 3); return (&x)[axis]; }
    __forceinline       T& operator []( const size_t axis )       { assert(axis < 3); return (&x)[axis]; }
  };

  ////////////////////////////////////////////////////////////////////////////////
  /// Unary Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline Vec3<T> operator +( const Vec3<T>& a ) { return Vec3<T>(+a.x, +a.y, +a.z); }
  template<typename T> __forceinline Vec3<T> operator -( const Vec3<T>& a ) { return Vec3<T>(-a.x, -a.y, -a.z); }
  template<typename T> __forceinline Vec3<T> abs       ( const Vec3<T>& a ) { return Vec3<T>(abs  (a.x), abs  (a.y), abs  (a.z)); }
  template<typename T> __forceinline Vec3<T> rcp       ( const Vec3<T>& a ) { return Vec3<T>(rcp  (a.x), rcp  (a.y), rcp  (a.z)); }
  template<typename T> __forceinline Vec3<T> rsqrt     ( const Vec3<T>& a ) { return Vec3<T>(rsqrt(a.x), rsqrt(a.y), rsqrt(a.z)); }
  template<typename T> __forceinline Vec3<T> sqrt      ( const Vec3<T>& a ) { return Vec3<T>(sqrt (a.x), sqrt (a.y), sqrt (a.z)); }
  template<typename T> __forceinline Vec3<T> zero_fix( const Vec3<T>& a ) {
    return Vec3<T>(select(a.x==0.0f,T(1E-10f),a.x),select(a.y==0.0f,T(1E-10f),a.y),select(a.z==0.0f,T(1E-10f),a.z));
  }
  template<typename T> __forceinline const Vec3<T> rcp_safe(const Vec3<T>& a) { return rcp(zero_fix(a)); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline Vec3<T> operator +( const Vec3<T>& a, const Vec3<T>& b ) { return Vec3<T>(a.x + b.x, a.y + b.y, a.z + b.z); }
  template<typename T> __forceinline Vec3<T> operator -( const Vec3<T>& a, const Vec3<T>& b ) { return Vec3<T>(a.x - b.x, a.y - b.y, a.z - b.z); }
  template<typename T> __forceinline Vec3<T> operator *( const Vec3<T>& a, const Vec3<T>& b ) { return Vec3<T>(a.x * b.x, a.y * b.y, a.z * b.z); }
  template<typename T> __forceinline Vec3<T> operator *( const       T& a, const Vec3<T>& b ) { return Vec3<T>(a   * b.x, a   * b.y, a   * b.z); }
  template<typename T> __forceinline Vec3<T> operator *( const Vec3<T>& a, const       T& b ) { return Vec3<T>(a.x * b  , a.y * b  , a.z * b  ); }
  template<typename T> __forceinline Vec3<T> operator /( const Vec3<T>& a, const       T& b ) { return Vec3<T>(a.x / b  , a.y / b  , a.z / b  ); }
  template<typename T> __forceinline Vec3<T> operator /( const       T& a, const Vec3<T>& b ) { return Vec3<T>(a   / b.x, a   / b.y, a   / b.z); }
  template<typename T> __forceinline Vec3<T> operator /( const Vec3<T>& a, const Vec3<T>& b ) { return Vec3<T>(a.x / b.x, a.y / b.y, a.z / b.z); }

  template<typename T> __forceinline Vec3<T> min(const Vec3<T>& a, const Vec3<T>& b) { return Vec3<T>(min(a.x, b.x), min(a.y, b.y), min(a.z, b.z)); }
  template<typename T> __forceinline Vec3<T> max(const Vec3<T>& a, const Vec3<T>& b) { return Vec3<T>(max(a.x, b.x), max(a.y, b.y), max(a.z, b.z)); }

  template<typename T> __forceinline Vec3<T> operator >>( const Vec3<T>& a, const int b ) { return Vec3<T>(a.x >> b, a.y >> b, a.z >> b); }
  template<typename T> __forceinline Vec3<T> operator <<( const Vec3<T>& a, const int b ) { return Vec3<T>(a.x << b, a.y << b, a.z << b); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Ternary Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline const Vec3<T> madd  ( const Vec3<T>& a, const Vec3<T>& b, const Vec3<T>& c) { return a*b+c; }
  template<typename T> __forceinline const Vec3<T> msub  ( const Vec3<T>& a, const Vec3<T>& b, const Vec3<T>& c) { return a*b-c; }
  template<typename T> __forceinline const Vec3<T> nmadd ( const Vec3<T>& a, const Vec3<T>& b, const Vec3<T>& c) { return -a*b-c;}
  template<typename T> __forceinline const Vec3<T> nmsub ( const Vec3<T>& a, const Vec3<T>& b, const Vec3<T>& c) { return c-a*b; }

  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline const Vec3<T>& operator +=( Vec3<T>& a, const T        b ) { a.x += b;   a.y += b;   a.z += b;   return a; }
  template<typename T> __forceinline const Vec3<T>& operator +=( Vec3<T>& a, const Vec3<T>& b ) { a.x += b.x; a.y += b.y; a.z += b.z; return a; }
  template<typename T> __forceinline const Vec3<T>& operator -=( Vec3<T>& a, const Vec3<T>& b ) { a.x -= b.x; a.y -= b.y; a.z -= b.z; return a; }
  template<typename T> __forceinline const Vec3<T>& operator *=( Vec3<T>& a, const       T& b ) { a.x *= b  ; a.y *= b  ; a.z *= b  ; return a; }
  template<typename T> __forceinline const Vec3<T>& operator /=( Vec3<T>& a, const       T& b ) { a.x /= b  ; a.y /= b  ; a.z /= b  ; return a; }

  ////////////////////////////////////////////////////////////////////////////////
  /// Reduction Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline T reduce_add( const Vec3<T>& a ) { return a.x + a.y + a.z; }
  template<typename T> __forceinline T reduce_mul( const Vec3<T>& a ) { return a.x * a.y * a.z; }
  template<typename T> __forceinline T reduce_min( const Vec3<T>& a ) { return min(a.x, a.y, a.z); }
  template<typename T> __forceinline T reduce_max( const Vec3<T>& a ) { return max(a.x, a.y, a.z); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline bool operator ==( const Vec3<T>& a, const Vec3<T>& b ) { return a.x == b.x && a.y == b.y && a.z == b.z; }
  template<typename T> __forceinline bool operator !=( const Vec3<T>& a, const Vec3<T>& b ) { return a.x != b.x || a.y != b.y || a.z != b.z; }
  template<typename T> __forceinline bool operator < ( const Vec3<T>& a, const Vec3<T>& b ) {
    if (a.x != b.x) return a.x < b.x;
    if (a.y != b.y) return a.y < b.y;
    if (a.z != b.z) return a.z < b.z;
    return false;
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Euclidian Space Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline T       dot      ( const Vec3<T>& a, const Vec3<T>& b ) { return madd(a.x,b.x,madd(a.y,b.y,a.z*b.z)); }
  template<typename T> __forceinline T       length   ( const Vec3<T>& a )                   { return sqrt(dot(a,a)); }
  template<typename T> __forceinline Vec3<T> normalize( const Vec3<T>& a )                   { return a*rsqrt(dot(a,a)); }

  template<typename T> __forceinline Vec3<T> normalize_safe( const Vec3<T>& a ) { 
    const T d = dot(a,a);
    return select(d == T( zero ), a ,  a*rsqrt(d) );
  }

  template<typename T> __forceinline T       distance ( const Vec3<T>& a, const Vec3<T>& b ) { return length(a-b); }
  template<typename T> __forceinline Vec3<T> cross    ( const Vec3<T>& a, const Vec3<T>& b ) { return Vec3<T>(msub(a.y,b.z,a.z*b.y), msub(a.z,b.x,a.x*b.z), msub(a.x,b.y,a.y*b.x)); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Select
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline Vec3<T> select ( bool s, const Vec3<T>& t, const Vec3<T>& f ) {
    return Vec3<T>(select(s,t.x,f.x),select(s,t.y,f.y),select(s,t.z,f.z));
  }

  template<typename T> __forceinline Vec3<T> select ( const Vec3<bool>& s, const Vec3<T>& t, const Vec3<T>& f ) {
    return Vec3<T>(select(s.x,t.x,f.x),select(s.y,t.y,f.y),select(s.z,t.z,f.z));
  }

  template<typename T> __forceinline Vec3<T> select ( const typename T::Mask& s, const Vec3<T>& t, const Vec3<T>& f ) {
    return Vec3<T>(select(s,t.x,f.x),select(s,t.y,f.y),select(s,t.z,f.z));
  }

  template<typename T> __forceinline int maxDim ( const Vec3<T>& a ) 
  { 
    if (a.x > a.y) {
      if (a.x > a.z) return 0; else return 2;
    } else {
      if (a.y > a.z) return 1; else return 2;
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators
  ////////////////////////////////////////////////////////////////////////////////
  template<typename T> __forceinline Vec3<bool> eq_mask( const Vec3<T>& a, const Vec3<T>& b ) { return Vec3<bool>(a.x==b.x,a.y==b.y,a.z==b.z); }
  template<typename T> __forceinline Vec3<bool> neq_mask(const Vec3<T>& a, const Vec3<T>& b ) { return Vec3<bool>(a.x!=b.x,a.y!=b.y,a.z!=b.z); }
  template<typename T> __forceinline Vec3<bool> lt_mask( const Vec3<T>& a, const Vec3<T>& b ) { return Vec3<bool>(a.x< b.x,a.y< b.y,a.z< b.z); }
  template<typename T> __forceinline Vec3<bool> le_mask( const Vec3<T>& a, const Vec3<T>& b ) { return Vec3<bool>(a.x<=b.x,a.y<=b.y,a.z<=b.z); }
  template<typename T> __forceinline Vec3<bool> gt_mask( const Vec3<T>& a, const Vec3<T>& b ) { return Vec3<bool>(a.x> b.x,a.y> b.y,a.z> b.z); }
  template<typename T> __forceinline Vec3<bool> ge_mask( const Vec3<T>& a, const Vec3<T>& b ) { return Vec3<bool>(a.x>=b.x,a.y>=b.y,a.z>=b.z); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Output Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> inline std::ostream& operator<<(std::ostream& cout, const Vec3<T>& a) {
    return cout << "(" << a.x << ", " << a.y << ", " << a.z << ")";
  }

  typedef Vec3<bool > Vec3b;
  typedef Vec3<int  > Vec3i;
  typedef Vec3<float> Vec3f;
}

#if defined(__MIC__)
#include "vec3ba_mic.h"
#include "vec3ia_mic.h"
#include "vec3fa_mic.h"
#else
#include "vec3ba.h" 
#include "vec3ia.h" 
#include "vec3fa.h" 
#endif

namespace embree 
{ 
  template<> __forceinline Vec3<float>::Vec3( const Vec3fa& a ) { x = a.x; y = a.y; z = a.z; }

#if defined (__SSE__)
  template<> __forceinline Vec3<ssef>::Vec3( const Vec3fa& a ) { 
    const ssef v = ssef(a); x = shuffle<0,0,0,0>(v); y = shuffle<1,1,1,1>(v); z = shuffle<2,2,2,2>(v); 
  }
  __forceinline Vec3<ssef> broadcast4f( const Vec3<ssef>& a, const size_t k ) {  
    return Vec3<ssef>(ssef::broadcast(&a.x[k]), ssef::broadcast(&a.y[k]), ssef::broadcast(&a.z[k]));
  }

  template<size_t i0, size_t i1, size_t i2, size_t i3> __forceinline const Vec3<ssef> shuffle( const Vec3<ssef>& b ) {
    return Vec3<ssef>(shuffle<i0,i1,i2,i3>(b.x),shuffle<i0,i1,i2,i3>(b.y),shuffle<i0,i1,i2,i3>(b.z));
  }

#endif

#if defined(__AVX__)
  template<> __forceinline Vec3<avxf>::Vec3( const Vec3fa& a ) {  
    x = a.x; y = a.y; z = a.z; 
  }
  __forceinline Vec3<ssef> broadcast4f( const Vec3<avxf>& a, const size_t k ) {  
    return Vec3<ssef>(ssef::broadcast(&a.x[k]), ssef::broadcast(&a.y[k]), ssef::broadcast(&a.z[k]));
  }
  __forceinline Vec3<avxf> broadcast8f( const Vec3<ssef>& a, const size_t k ) {  
    return Vec3<avxf>(avxf::broadcast(&a.x[k]), avxf::broadcast(&a.y[k]), avxf::broadcast(&a.z[k]));
  }
  __forceinline Vec3<avxf> broadcast8f( const Vec3<avxf>& a, const size_t k ) {  
    return Vec3<avxf>(avxf::broadcast(&a.x[k]), avxf::broadcast(&a.y[k]), avxf::broadcast(&a.z[k]));
  }

  template<size_t i0, size_t i1, size_t i2, size_t i3> __forceinline const Vec3<avxf> shuffle( const Vec3<avxf>& b ) {
    return Vec3<avxf>(shuffle<i0,i1,i2,i3>(b.x),shuffle<i0,i1,i2,i3>(b.y),shuffle<i0,i1,i2,i3>(b.z));
  }

#endif

#if defined(__MIC__)
  //template<> __forceinline Vec3<ssef>::Vec3( const Vec3fa& a ) : x(a.x), y(a.y), z(a.z) {}
  template<> __forceinline Vec3<mic_f>::Vec3( const Vec3fa& a ) : x(a.x), y(a.y), z(a.z) {}
#endif
}

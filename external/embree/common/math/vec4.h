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

namespace embree
{
  ////////////////////////////////////////////////////////////////////////////////
  /// Generic 4D vector Class
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> struct Vec4
  {
    T x, y, z, w;

    typedef T Scalar;
    enum { N = 4 };

    ////////////////////////////////////////////////////////////////////////////////
    /// Construction
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vec4    ( )                  { }
    __forceinline Vec4    ( const Vec4& other ) { x = other.x; y = other.y; z = other.z; w = other.w; }
    __forceinline Vec4    ( const Vec3fa& other );

    template<typename T1> __forceinline Vec4( const Vec4<T1>& a ) : x(T(a.x)), y(T(a.y)), z(T(a.z)), w(T(a.w)) {}
    template<typename T1> __forceinline Vec4& operator =(const Vec4<T1>& other) { x = other.x; y = other.y; z = other.z; w = other.w; return *this; }

    __forceinline explicit Vec4( const T& a                                     ) : x(a), y(a), z(a), w(a) {}
    __forceinline explicit Vec4( const T& x, const T& y, const T& z, const T& w ) : x(x), y(y), z(z), w(w) {}
    __forceinline explicit Vec4( const T* const a, const size_t stride = 1     ) : x(a[0]), y(a[stride]), z(a[2*stride]), w(a[3*stride]) {}

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vec4( ZeroTy   ) : x(zero), y(zero), z(zero), w(zero) {}
    __forceinline Vec4( OneTy    ) : x(one),  y(one),  z(one),  w(one) {}
    __forceinline Vec4( PosInfTy ) : x(pos_inf), y(pos_inf), z(pos_inf), w(pos_inf) {}
    __forceinline Vec4( NegInfTy ) : x(neg_inf), y(neg_inf), z(neg_inf), w(neg_inf) {}

    __forceinline const T& operator []( const size_t axis ) const { assert(axis < 4); return (&x)[axis]; }
    __forceinline       T& operator []( const size_t axis )       { assert(axis < 4); return (&x)[axis]; }
  };

  ////////////////////////////////////////////////////////////////////////////////
  /// Unary Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline Vec4<T> operator +( const Vec4<T>& a ) { return Vec4<T>(+a.x, +a.y, +a.z, +a.w); }
  template<typename T> __forceinline Vec4<T> operator -( const Vec4<T>& a ) { return Vec4<T>(-a.x, -a.y, -a.z, -a.w); }
  template<typename T> __forceinline Vec4<T> abs       ( const Vec4<T>& a ) { return Vec4<T>(abs  (a.x), abs  (a.y), abs  (a.z), abs  (a.w)); }
  template<typename T> __forceinline Vec4<T> rcp       ( const Vec4<T>& a ) { return Vec4<T>(rcp  (a.x), rcp  (a.y), rcp  (a.z), rcp  (a.w)); }
  template<typename T> __forceinline Vec4<T> rsqrt     ( const Vec4<T>& a ) { return Vec4<T>(rsqrt(a.x), rsqrt(a.y), rsqrt(a.z), rsqrt(a.w)); }
  template<typename T> __forceinline Vec4<T> sqrt      ( const Vec4<T>& a ) { return Vec4<T>(sqrt (a.x), sqrt (a.y), sqrt (a.z), sqrt (a.w)); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline Vec4<T> operator +( const Vec4<T>& a, const Vec4<T>& b ) { return Vec4<T>(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w); }
  template<typename T> __forceinline Vec4<T> operator -( const Vec4<T>& a, const Vec4<T>& b ) { return Vec4<T>(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w); }
  template<typename T> __forceinline Vec4<T> operator *( const Vec4<T>& a, const Vec4<T>& b ) { return Vec4<T>(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w); }
  template<typename T> __forceinline Vec4<T> operator *( const       T& a, const Vec4<T>& b ) { return Vec4<T>(a   * b.x, a   * b.y, a   * b.z, a   * b.w); }
  template<typename T> __forceinline Vec4<T> operator *( const Vec4<T>& a, const       T& b ) { return Vec4<T>(a.x * b  , a.y * b  , a.z * b  , a.w * b  ); }
  template<typename T> __forceinline Vec4<T> operator /( const Vec4<T>& a, const Vec4<T>& b ) { return Vec4<T>(a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w); }
  template<typename T> __forceinline Vec4<T> operator /( const Vec4<T>& a, const       T& b ) { return Vec4<T>(a.x / b  , a.y / b  , a.z / b  , a.w / b  ); }
  template<typename T> __forceinline Vec4<T> operator /( const       T& a, const Vec4<T>& b ) { return Vec4<T>(a   / b.x, a   / b.y, a   / b.z, a   / b.w); }

  template<typename T> __forceinline Vec4<T> min(const Vec4<T>& a, const Vec4<T>& b) { return Vec4<T>(min(a.x, b.x), min(a.y, b.y), min(a.z, b.z), min(a.w, b.w)); }
  template<typename T> __forceinline Vec4<T> max(const Vec4<T>& a, const Vec4<T>& b) { return Vec4<T>(max(a.x, b.x), max(a.y, b.y), max(a.z, b.z), max(a.w, b.w)); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline Vec4<T>& operator +=( Vec4<T>& a, const Vec4<T>& b ) { a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w; return a; }
  template<typename T> __forceinline Vec4<T>& operator -=( Vec4<T>& a, const Vec4<T>& b ) { a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w; return a; }
  template<typename T> __forceinline Vec4<T>& operator *=( Vec4<T>& a, const       T& b ) { a.x *= b  ; a.y *= b  ; a.z *= b  ; a.w *= b  ; return a; }
  template<typename T> __forceinline Vec4<T>& operator /=( Vec4<T>& a, const       T& b ) { a.x /= b  ; a.y /= b  ; a.z /= b  ; a.w /= b  ; return a; }

  ////////////////////////////////////////////////////////////////////////////////
  /// Reduction Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline T reduce_add( const Vec4<T>& a ) { return a.x + a.y + a.z + a.w; }
  template<typename T> __forceinline T reduce_mul( const Vec4<T>& a ) { return a.x * a.y * a.z * a.w; }
  template<typename T> __forceinline T reduce_min( const Vec4<T>& a ) { return min(a.x, a.y, a.z, a.w); }
  template<typename T> __forceinline T reduce_max( const Vec4<T>& a ) { return max(a.x, a.y, a.z, a.w); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline bool operator ==( const Vec4<T>& a, const Vec4<T>& b ) { return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w; }
  template<typename T> __forceinline bool operator !=( const Vec4<T>& a, const Vec4<T>& b ) { return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w; }
  template<typename T> __forceinline bool operator < ( const Vec4<T>& a, const Vec4<T>& b ) {
    if (a.x != b.x) return a.x < b.x;
    if (a.y != b.y) return a.y < b.y;
    if (a.z != b.z) return a.z < b.z;
    if (a.w != b.w) return a.w < b.w;
    return false;
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Euclidian Space Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline T       dot      ( const Vec4<T>& a, const Vec4<T>& b ) { return a.x*b.x + a.y*b.y + a.z*b.z + a.w*b.w; }
  template<typename T> __forceinline T       length   ( const Vec4<T>& a )                   { return sqrt(dot(a,a)); }
  template<typename T> __forceinline Vec4<T> normalize( const Vec4<T>& a )                   { return a*rsqrt(dot(a,a)); }
  template<typename T> __forceinline T       distance ( const Vec4<T>& a, const Vec4<T>& b ) { return length(a-b); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Select
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline Vec4<T> select ( bool s, const Vec4<T>& t, const Vec4<T>& f ) {
    return Vec4<T>(select(s,t.x,f.x),select(s,t.y,f.y),select(s,t.z,f.z),select(s,t.w,f.w));
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Output Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> inline std::ostream& operator<<(std::ostream& cout, const Vec4<T>& a) {
    return cout << "(" << a.x << ", " << a.y << ", " << a.z << ", " << a.w << ")";
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Default template instantiations
  ////////////////////////////////////////////////////////////////////////////////

  typedef Vec4<bool         > Vec4b;
  typedef Vec4<unsigned char> Vec4uc;
  typedef Vec4<int          > Vec4i;
  typedef Vec4<float        > Vec4f;
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
  template<> __forceinline Vec4<float>::Vec4( const Vec3fa& a ) { x = a.x; y = a.y; z = a.z; w = a.w; }

#if defined (__SSE__) 
  template<> __forceinline Vec4<ssef>::Vec4( const Vec3fa& a ) { 
    const ssef v = ssef(a); x = shuffle<0,0,0,0>(v); y = shuffle<1,1,1,1>(v); z = shuffle<2,2,2,2>(v); w = shuffle<3,3,3,3>(v); 
  }
  __forceinline Vec4<ssef> broadcast4f( const Vec4<ssef>& a, const size_t k ) {  
    return Vec4<ssef>(ssef::broadcast(&a.x[k]), ssef::broadcast(&a.y[k]), ssef::broadcast(&a.z[k]), ssef::broadcast(&a.w[k]));
  }
#endif

#if defined(__AVX__)
  template<> __forceinline Vec4<avxf>::Vec4( const Vec3fa& a ) {  
    x = a.x; y = a.y; z = a.z; w = a.w; 
  }
  __forceinline Vec4<ssef> broadcast4f( const Vec4<avxf>& a, const size_t k ) {  
    return Vec4<ssef>(ssef::broadcast(&a.x[k]), ssef::broadcast(&a.y[k]), ssef::broadcast(&a.z[k]), ssef::broadcast(&a.w[k]));
  }
  __forceinline Vec4<avxf> broadcast8f( const Vec4<ssef>& a, const size_t k ) {  
    return Vec4<avxf>(avxf::broadcast(&a.x[k]), avxf::broadcast(&a.y[k]), avxf::broadcast(&a.z[k]), avxf::broadcast(&a.w[k]));
  }
  __forceinline Vec4<avxf> broadcast8f( const Vec4<avxf>& a, const size_t k ) {  
    return Vec4<avxf>(avxf::broadcast(&a.x[k]), avxf::broadcast(&a.y[k]), avxf::broadcast(&a.z[k]), avxf::broadcast(&a.w[k]));
  }
#endif

#if defined(__MIC__)
  //template<> __forceinline Vec4<ssef>::Vec4( const Vec3fa& a ) : x(a.x), y(a.y), z(a.z), w(a.w) {}
  template<> __forceinline Vec4<mic_f>::Vec4( const Vec3fa& a ) : x(a.x), y(a.y), z(a.z), w(a.w) {}
#endif
}

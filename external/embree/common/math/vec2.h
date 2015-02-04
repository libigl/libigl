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

namespace embree
{
  ////////////////////////////////////////////////////////////////////////////////
  /// Generic 2D vector Class
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> struct Vec2
  {
    T x, y;

    typedef T Scalar;
    enum { N = 2 };

    ////////////////////////////////////////////////////////////////////////////////
    /// Construction
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vec2     ( )                  { }
    __forceinline Vec2     ( const Vec2& other ) { x = other.x; y = other.y; }
    template<typename T1> __forceinline Vec2( const Vec2<T1>& a ) : x(T(a.x)), y(T(a.y)) {}
    template<typename T1> __forceinline Vec2& operator =( const Vec2<T1>& other ) { x = other.x; y = other.y; return *this; }

    __forceinline explicit Vec2( const T& a             ) : x(a), y(a) {}
    __forceinline explicit Vec2( const T& x, const T& y ) : x(x), y(y) {}
    __forceinline explicit Vec2( const T* const a, const ssize_t stride = 1 ) : x(a[0]), y(a[stride]) {}

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vec2( ZeroTy   ) : x(zero), y(zero) {}
    __forceinline Vec2( OneTy    ) : x(one),  y(one) {}
    __forceinline Vec2( PosInfTy ) : x(pos_inf), y(pos_inf) {}
    __forceinline Vec2( NegInfTy ) : x(neg_inf), y(neg_inf) {}

    __forceinline const T& operator []( const size_t axis ) const { assert(axis < 2); return (&x)[axis]; }
    __forceinline       T& operator []( const size_t axis )       { assert(axis < 2); return (&x)[axis]; }
  };

  ////////////////////////////////////////////////////////////////////////////////
  /// Unary Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline Vec2<T> operator +( const Vec2<T>& a ) { return Vec2<T>(+a.x, +a.y); }
  template<typename T> __forceinline Vec2<T> operator -( const Vec2<T>& a ) { return Vec2<T>(-a.x, -a.y); }
  template<typename T> __forceinline Vec2<T> abs       ( const Vec2<T>& a ) { return Vec2<T>(abs  (a.x), abs  (a.y)); }
  template<typename T> __forceinline Vec2<T> rcp       ( const Vec2<T>& a ) { return Vec2<T>(rcp  (a.x), rcp  (a.y)); }
  template<typename T> __forceinline Vec2<T> rsqrt     ( const Vec2<T>& a ) { return Vec2<T>(rsqrt(a.x), rsqrt(a.y)); }
  template<typename T> __forceinline Vec2<T> sqrt      ( const Vec2<T>& a ) { return Vec2<T>(sqrt (a.x), sqrt (a.y)); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline Vec2<T> operator +( const Vec2<T>& a, const Vec2<T>& b ) { return Vec2<T>(a.x + b.x, a.y + b.y); }
  template<typename T> __forceinline Vec2<T> operator -( const Vec2<T>& a, const Vec2<T>& b ) { return Vec2<T>(a.x - b.x, a.y - b.y); }
  template<typename T> __forceinline Vec2<T> operator *( const Vec2<T>& a, const Vec2<T>& b ) { return Vec2<T>(a.x * b.x, a.y * b.y); }
  template<typename T> __forceinline Vec2<T> operator *( const       T& a, const Vec2<T>& b ) { return Vec2<T>(a   * b.x, a   * b.y); }
  template<typename T> __forceinline Vec2<T> operator *( const Vec2<T>& a, const       T& b ) { return Vec2<T>(a.x * b  , a.y * b  ); }
  template<typename T> __forceinline Vec2<T> operator /( const Vec2<T>& a, const Vec2<T>& b ) { return Vec2<T>(a.x / b.x, a.y / b.y); }
  template<typename T> __forceinline Vec2<T> operator /( const Vec2<T>& a, const       T& b ) { return Vec2<T>(a.x / b  , a.y / b  ); }
  template<typename T> __forceinline Vec2<T> operator /( const       T& a, const Vec2<T>& b ) { return Vec2<T>(a   / b.x, a   / b.y); }

  template<typename T> __forceinline Vec2<T> min(const Vec2<T>& a, const Vec2<T>& b) { return Vec2<T>(min(a.x, b.x), min(a.y, b.y)); }
  template<typename T> __forceinline Vec2<T> max(const Vec2<T>& a, const Vec2<T>& b) { return Vec2<T>(max(a.x, b.x), max(a.y, b.y)); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline Vec2<T>& operator +=( Vec2<T>& a, const Vec2<T>& b ) { a.x += b.x; a.y += b.y; return a; }
  template<typename T> __forceinline Vec2<T>& operator -=( Vec2<T>& a, const Vec2<T>& b ) { a.x -= b.x; a.y -= b.y; return a; }
  template<typename T> __forceinline Vec2<T>& operator *=( Vec2<T>& a, const       T& b ) { a.x *= b  ; a.y *= b  ; return a; }
  template<typename T> __forceinline Vec2<T>& operator /=( Vec2<T>& a, const       T& b ) { a.x /= b  ; a.y /= b  ; return a; }

  ////////////////////////////////////////////////////////////////////////////////
  /// Reduction Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline T reduce_add( const Vec2<T>& a ) { return a.x + a.y; }
  template<typename T> __forceinline T reduce_mul( const Vec2<T>& a ) { return a.x * a.y; }
  template<typename T> __forceinline T reduce_min( const Vec2<T>& a ) { return min(a.x, a.y); }
  template<typename T> __forceinline T reduce_max( const Vec2<T>& a ) { return max(a.x, a.y); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline bool operator ==( const Vec2<T>& a, const Vec2<T>& b ) { return a.x == b.x && a.y == b.y; }
  template<typename T> __forceinline bool operator !=( const Vec2<T>& a, const Vec2<T>& b ) { return a.x != b.x || a.y != b.y; }
  template<typename T> __forceinline bool operator < ( const Vec2<T>& a, const Vec2<T>& b ) {
    if (a.x != b.x) return a.x < b.x;
    if (a.y != b.y) return a.y < b.y;
    return false;
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Euclidian Space Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline T       dot      ( const Vec2<T>& a, const Vec2<T>& b ) { return a.x*b.x + a.y*b.y; }
  template<typename T> __forceinline T       length   ( const Vec2<T>& a )                   { return sqrt(dot(a,a)); }
  template<typename T> __forceinline Vec2<T> normalize( const Vec2<T>& a )                   { return a*rsqrt(dot(a,a)); }
  template<typename T> __forceinline T       distance ( const Vec2<T>& a, const Vec2<T>& b ) { return length(a-b); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Select
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline Vec2<T> select ( bool s, const Vec2<T>& t, const Vec2<T>& f ) {
    return Vec2<T>(select(s,t.x,f.x),select(s,t.y,f.y));
  }

  template<typename T> __forceinline Vec2<T> select ( const typename T::Mask& s, const Vec2<T>& t, const Vec2<T>& f ) {
    return Vec2<T>(select(s,t.x,f.x),select(s,t.y,f.y));
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Output Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> inline std::ostream& operator<<(std::ostream& cout, const Vec2<T>& a) {
    return cout << "(" << a.x << ", " << a.y << ")";
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Default template instantiations
  ////////////////////////////////////////////////////////////////////////////////

  typedef Vec2<bool > Vec2b;
  typedef Vec2<int  > Vec2i;
  typedef Vec2<float> Vec2f;
}

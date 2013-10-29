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

#ifndef __EMBREE_VEC3_H__
#define __EMBREE_VEC3_H__

#include "sys/platform.h"
#include "math/math.h"

//!< use SSE for short vectors only in 64 bit mode (32 bit mode causes alignment issues)
#if defined(__MIC__) && 1
#  include "vector3b_mic.h" 
#  include "vector3i_mic.h" 
#  include "vector3f_mic.h" 
#  define EMBREE_FAST_VEC3
#elif defined(__X86_64__) && 1
#  include "vector3b_sse.h" 
#  include "vector3i_sse.h" 
#  include "vector3f_sse.h" 
#  define EMBREE_FAST_VEC3
#endif

namespace embree
{
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

    __forceinline Vec3     ( )                  { }
    __forceinline Vec3     ( const Vec3& other ) { x = other.x; y = other.y; z = other.z; }
    template<typename T1> __forceinline Vec3( const Vec3<T1>& a ) : x(T(a.x)), y(T(a.y)), z(T(a.z)) {}

#if defined EMBREE_FAST_VEC3
    __forceinline Vec3( const Vector3f& a ) : x(a.x), y(a.y), z(a.z) {}
    __forceinline Vec3( const Vec3fa& a ) : x(a.x), y(a.y), z(a.z) {}
#endif

    template<typename T1> __forceinline Vec3& operator =(const Vec3<T1>& other) { x = other.x; y = other.y; z = other.z; return *this; }

    __forceinline Vec3( const T& a                         ) : x(a), y(a), z(a) {}
    __forceinline explicit Vec3( const T& x, const T& y, const T& z ) : x(x), y(y), z(z) {}
    __forceinline explicit Vec3( const T* const a, const size_t stride = 1 ) : x(a[0]), y(a[stride]), z(a[2*stride]) {}

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

#if defined __EMBREE_VEC3F_SSE_H__
  template<> __forceinline Vec3<ssef>::Vec3( const Vector3f& a ) { const ssef v = ssef(a); x = shuffle<0,0,0,0>(v); y = shuffle<1,1,1,1>(v); z = shuffle<2,2,2,2>(v); }
  template<> __forceinline Vec3<ssef>::Vec3( const Vec3fa& a ) { const ssef v = ssef(Vector3f(a)); x = shuffle<0,0,0,0>(v); y = shuffle<1,1,1,1>(v); z = shuffle<2,2,2,2>(v); }
#endif

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
  /// Output Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> inline std::ostream& operator<<(std::ostream& cout, const Vec3<T>& a) {
    return cout << "(" << a.x << ", " << a.y << ", " << a.z << ")";
  }

  typedef Vec3<bool > Vec3b;
  typedef Vec3<int  > Vec3i;
  typedef Vec3<float> Vec3f;
}

////////////////////////////////////////////////////////////////////////////////
/// Comparison Operators
////////////////////////////////////////////////////////////////////////////////
namespace embree {
  template<typename T> __forceinline Vec3<bool> eq_mask( const Vec3<T>& a, const Vec3<T>& b ) { return Vec3<bool>(a.x==b.x,a.y==b.y,a.z==b.z); }
  template<typename T> __forceinline Vec3<bool> neq_mask(const Vec3<T>& a, const Vec3<T>& b ) { return Vec3<bool>(a.x!=b.x,a.y!=b.y,a.z!=b.z); }
  template<typename T> __forceinline Vec3<bool> lt_mask( const Vec3<T>& a, const Vec3<T>& b ) { return Vec3<bool>(a.x< b.x,a.y< b.y,a.z< b.z); }
  template<typename T> __forceinline Vec3<bool> le_mask( const Vec3<T>& a, const Vec3<T>& b ) { return Vec3<bool>(a.x<=b.x,a.y<=b.y,a.z<=b.z); }
  template<typename T> __forceinline Vec3<bool> gt_mask( const Vec3<T>& a, const Vec3<T>& b ) { return Vec3<bool>(a.x> b.x,a.y> b.y,a.z> b.z); }
  template<typename T> __forceinline Vec3<bool> ge_mask( const Vec3<T>& a, const Vec3<T>& b ) { return Vec3<bool>(a.x>=b.x,a.y>=b.y,a.z>=b.z); }
}

#if !defined(EMBREE_FAST_VEC3)
namespace embree 
{
  class Vec3fa : public Vec3f { 
  public:
    __forceinline Vec3fa () {}
    __forceinline Vec3fa (const float v) : Vec3f(v), a(0) {}
    __forceinline Vec3fa (const float x, const float y, const float z) : Vec3f(x,y,z), a(0) {}
    __forceinline Vec3fa (const Vec3f& v) : Vec3f(v), a(0) {}
  public:
    int a; 
  };

  typedef Vec3<bool > Vector3b;
  typedef Vec3<int  > Vector3i;
  typedef Vec3<float> Vector3f;
}
#endif

#endif

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

#ifndef __EMBREE_COL3_H__
#define __EMBREE_COL3_H__

#include "sys/platform.h"
#include "math/math.h"

namespace embree
{
  ////////////////////////////////////////////////////////////////////////////////
  /// RGB Color Class
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> struct Col3
  {
    T r, g, b;

    ////////////////////////////////////////////////////////////////////////////////
    /// Construction
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Col3           ( )                   { }
    __forceinline Col3           ( const Col3& other ) { r = other.r; g = other.g; b = other.b; }
    __forceinline Col3& operator=( const Col3& other ) { r = other.r; g = other.g; b = other.b; return *this; }
    template<typename T1> __forceinline explicit Col3( const Col3<T1>& a   ) : r(T(a.r)), g(T(a.g)), b(T(a.b)) {}

    __forceinline explicit Col3 (const T& v)                         : r(v), g(v), b(v) {}
    __forceinline explicit Col3 (const T* v, int stride = 1)         : r(v[0*stride]), g(v[1*stride]), b(v[2*stride]) {}
    __forceinline          Col3 (const T& r, const T& g, const T& b) : r(r), g(g), b(b) {}

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Col3 (ZeroTy)   : r(zero)   , g(zero)   , b(zero)    {}
    __forceinline Col3 (OneTy)    : r(one)    , g(one)    , b(one)     {}
    __forceinline Col3 (PosInfTy) : r(pos_inf), g(pos_inf), b(pos_inf) {}
    __forceinline Col3 (NegInfTy) : r(neg_inf), g(neg_inf), b(neg_inf) {}
  };

  ////////////////////////////////////////////////////////////////////////////////
  /// Unary Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline const Col3<T> operator+ (const Col3<T>& v) { return Col3<T>(+v.r,+v.g,+v.b); }
  template<typename T> __forceinline const Col3<T> operator- (const Col3<T>& v) { return Col3<T>(-v.r,-v.g,-v.b); }
  template<typename T> __forceinline const Col3<T> abs       (const Col3<T>& a) { return Col3<T>(abs(a.r),abs(a.g),abs(a.b)); }
  template<typename T> __forceinline const Col3<T> rcp       (const Col3<T>& a) { return Col3<T>(rcp(a.r),rcp(a.g),rcp(a.b)); }
  template<typename T> __forceinline const Col3<T> rsqrt     (const Col3<T>& a) { return Col3<T>(rsqrt(a.r),rsqrt(a.g),rsqrt(a.b)); }
  template<typename T> __forceinline const Col3<T> sqrt      (const Col3<T>& a) { return Col3<T>(sqrt(a.r),sqrt(a.g),sqrt(a.b)); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline const Col3<T> operator+(const Col3<T>& a, const Col3<T>& b) { return Col3<T>(a.r+b.r,a.g+b.g,a.b+b.b); }
  template<typename T> __forceinline const Col3<T> operator-(const Col3<T>& a, const Col3<T>& b) { return Col3<T>(a.r-b.r,a.g-b.g,a.b-b.b); }
  template<typename T> __forceinline const Col3<T> operator*(const T& a,       const Col3<T>& b) { return Col3<T>(a*b.r,a*b.g,a*b.b); }
  template<typename T> __forceinline const Col3<T> operator*(const Col3<T>& a, const T& b      ) { return Col3<T>(a.r*b,a.g*b,a.b*b); }
  template<typename T> __forceinline const Col3<T> operator*(const Col3<T>& a, const Col3<T>& b) { return Col3<T>(a.r*b.r,a.g*b.g,a.b*b.b); }
  template<typename T> __forceinline const Col3<T> operator/(const Col3<T>& a, const Col3<T>& b) { return Col3<T>(a.r/b.r,a.g/b.g,a.b/b.b); }
  template<typename T> __forceinline const Col3<T> operator/(const Col3<T>& a, const T& b      ) { return Col3<T>(a.r/b,a.g/b,a.b/b); }

  template<typename T> __forceinline const Col3<T> min(const Col3<T> a, const Col3<T> b) { return Col3<T>(min(a.r,b.r),min(a.g,b.g),min(a.b,b.b)); }
  template<typename T> __forceinline const Col3<T> max(const Col3<T> a, const Col3<T> b) { return Col3<T>(max(a.r,b.r),max(a.g,b.g),max(a.b,b.b)); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline const Col3<T> operator+=(Col3<T>& a, const Col3<T>& b) { return a = a + b; }
  template<typename T> __forceinline const Col3<T> operator-=(Col3<T>& a, const Col3<T>& b) { return a = a - b; }
  template<typename T> __forceinline const Col3<T> operator*=(Col3<T>& a, const Col3<T>& b) { return a = a * b; }
  template<typename T> __forceinline const Col3<T> operator/=(Col3<T>& a, const Col3<T>& b) { return a = a / b; }
  template<typename T> __forceinline const Col3<T> operator*=(Col3<T>& a, const T& b)       { return a = a * b; }
  template<typename T> __forceinline const Col3<T> operator/=(Col3<T>& a, const T& b)       { return a = a / b; }

  ////////////////////////////////////////////////////////////////////////////////
  /// Reduction Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline const T reduce_add(const Col3<T>& a) { return a.r+a.g+a.b; }
  template<typename T> __forceinline const T reduce_mul(const Col3<T>& a) { return a.r*a.g*a.b; }
  template<typename T> __forceinline const T reduce_min(const Col3<T>& a) { return min(a.r,a.g,a.b); }
  template<typename T> __forceinline const T reduce_max(const Col3<T>& a) { return max(a.r,a.g,a.b); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline bool operator ==(const Col3<T>& a, const Col3<T>& b) { return a.r == b.r && a.g == b.g && a.b == b.b; }
  template<typename T> __forceinline bool operator !=(const Col3<T>& a, const Col3<T>& b) { return a.r != b.r || a.g != b.g || a.b != b.b; }
  template<typename T> __forceinline bool operator < (const Col3<T>& a, const Col3<T>& b ) {
    if (a.r != b.r) return a.r < b.r;
    if (a.g != b.g) return a.g < b.g;
    if (a.b != b.b) return a.b < b.b;
    return false;
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Select
  ////////////////////////////////////////////////////////////////////////////////

  template<typename T> __forceinline Col3<T> select ( bool s, const Col3<T>& t, const Col3<T>& f )
  { return Col3<T>(select(s,t.r,f.r),select(s,t.g,f.g),select(s,t.b,f.b)); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Special Operators
  ////////////////////////////////////////////////////////////////////////////////

  /*! computes luminance of a color */
  template<typename T> __forceinline const T luminance (const Col3<T>& a) { return 0.212671f*a.r + 0.715160f*a.g + 0.072169f*a.b; }

  template<typename T> __forceinline Col3<T> exp (const Col3<T>& a) { return Col3<T>(exp(a.r),exp(a.g),exp(a.b)); }
  template<typename T> __forceinline Col3<T> log (const Col3<T>& a) { return Col3<T>(log(a.r),log(a.g),log(a.b)); }
  template<typename T> __forceinline Col3<T> pow (const Col3<T>& a, float e) { return exp(log(max(Col3<T>(1E-10f),a))*e); }

  /*! output operator */
  template<typename T> inline std::ostream& operator<<(std::ostream& cout, const Col3<T>& a) {
    return cout << "(" << a.r << ", " << a.g << ", " << a.b << ")";
  }

  /*! default template instantiations */
  typedef Col3<unsigned char> Col3c;
  typedef Col3<float> Col3f;
}

#if defined(__X86_64__)     //!< 32 bit mode causes aligmnent issues
#  include "math/col3f_sse.h" //!< comment out for not using SSE colors
#endif

#endif

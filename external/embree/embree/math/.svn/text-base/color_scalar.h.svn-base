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

#ifndef __EMBREE_COLOR_SCALAR_H__
#define __EMBREE_COLOR_SCALAR_H__

#include "sys/platform.h"
#include "math/math.h"

#include "col3.h"
#include "col4.h"

namespace embree
{
  ////////////////////////////////////////////////////////////////////////////////
  /// RGBA Color Class
  ////////////////////////////////////////////////////////////////////////////////
  
  struct Color4
  {
    float r, g, b, a;

    ////////////////////////////////////////////////////////////////////////////////
    /// Construction
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Color4 () {}

    __forceinline explicit Color4 (const float& v) : r(v), g(v), b(v), a(v) {}
    __forceinline          Color4 (const float& r, const float& g, const float& b, const float& a) : r(r), g(g), b(b), a(a) {}

    __forceinline explicit Color4 ( const Col3c& other ) { r = other.r*one_over_255; g = other.g*one_over_255; b = other.b*one_over_255; a = 1.0f; }
    __forceinline explicit Color4 ( const Col3f& other ) { r = other.r; g = other.g; b = other.b; a = 1.0f; }
    __forceinline explicit Color4 ( const Col4c& other ) { r = other.r*one_over_255; g = other.g*one_over_255; b = other.b*one_over_255; a = other.a*one_over_255; }
    __forceinline explicit Color4 ( const Col4f& other ) { r = other.r; g = other.g; b = other.b; a = other.a; }
    
    __forceinline Color4           ( const Color4& other ) { r = other.r; g = other.g; b = other.b; a = other.a; }
    __forceinline Color4& operator=( const Color4& other ) { r = other.r; g = other.g; b = other.b; a = other.a; return *this; }

    ////////////////////////////////////////////////////////////////////////////////
    /// Set
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline void set(Col3f& d) const { d.r = r; d.g = g; d.b = b; }
    __forceinline void set(Col4f& d) const { d.r = r; d.g = g; d.b = b; d.a = a; }
    __forceinline void set(Col3c& d) const { d.r = char(clamp(r)*255.0f); d.g = char(clamp(g)*255.0f); d.b = char(clamp(b)*255.0f); }
    __forceinline void set(Col4c& d) const { d.r = char(clamp(r)*255.0f); d.g = char(clamp(g)*255.0f); d.b = char(clamp(b)*255.0f); d.a = char(clamp(a)*255.0f); }

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Color4 (ZeroTy)   : r(zero)   , g(zero)   , b(zero)   , a(zero)    {}
    __forceinline Color4 (OneTy)    : r(one)    , g(one)    , b(one)    , a(one)     {}
    __forceinline Color4 (PosInfTy) : r(pos_inf), g(pos_inf), b(pos_inf), a(pos_inf) {}
    __forceinline Color4 (NegInfTy) : r(neg_inf), g(neg_inf), b(neg_inf), a(neg_inf) {}
  };

  ////////////////////////////////////////////////////////////////////////////////
  /// RGB Color Class
  ////////////////////////////////////////////////////////////////////////////////

  struct Color
  {
    float r, g, b;

    ////////////////////////////////////////////////////////////////////////////////
    /// Construction
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Color () {}

    __forceinline explicit Color (const float& v) : r(v), g(v), b(v) {}
    __forceinline          Color (const float& r, const float& g, const float& b) : r(r), g(g), b(b) {}

    __forceinline Color           ( const Color& other ) { r = other.r; g = other.g; b = other.b; }
    __forceinline Color& operator=( const Color& other ) { r = other.r; g = other.g; b = other.b; return *this; }

    __forceinline Color           ( const Color4& other ) { r = other.r; g = other.g; b = other.b; }
    __forceinline Color& operator=( const Color4& other ) { r = other.r; g = other.g; b = other.b; return *this; }

    ////////////////////////////////////////////////////////////////////////////////
    /// Set
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline void set(Col3f& d) const { d.r = r; d.g = g; d.b = b; }
    __forceinline void set(Col4f& d) const { d.r = r; d.g = g; d.b = b; d.a = 1.0f; }
    __forceinline void set(Col3c& d) const { d.r = char(clamp(r)*255.0f); d.g = char(clamp(g)*255.0f); d.b = char(clamp(b)*255.0f); }
    __forceinline void set(Col4c& d) const { d.r = char(clamp(r)*255.0f); d.g = char(clamp(g)*255.0f); d.b = char(clamp(b)*255.0f); d.a = 1.0f; }

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Color (ZeroTy)   : r(zero)   , g(zero)   , b(zero)    {}
    __forceinline Color (OneTy)    : r(one)    , g(one)    , b(one)     {}
    __forceinline Color (PosInfTy) : r(pos_inf), g(pos_inf), b(pos_inf) {}
    __forceinline Color (NegInfTy) : r(neg_inf), g(neg_inf), b(neg_inf) {}
  };

  ////////////////////////////////////////////////////////////////////////////////
  /// Unary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Color operator+ (const Color& v) { return Color(+v.r,+v.g,+v.b); }
  __forceinline const Color operator- (const Color& v) { return Color(-v.r,-v.g,-v.b); }
  __forceinline const Color abs       (const Color& a) { return Color(abs(a.r),abs(a.g),abs(a.b)); }
  __forceinline const Color rcp       (const Color& a) { return Color(rcp(a.r),rcp(a.g),rcp(a.b)); }
  __forceinline const Color rsqrt     (const Color& a) { return Color(rsqrt(a.r),rsqrt(a.g),rsqrt(a.b)); }
  __forceinline const Color sqrt      (const Color& a) { return Color(sqrt(a.r),sqrt(a.g),sqrt(a.b)); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Color operator+(const Color& a, const Color& b) { return Color(a.r+b.r,a.g+b.g,a.b+b.b); }
  __forceinline const Color operator-(const Color& a, const Color& b) { return Color(a.r-b.r,a.g-b.g,a.b-b.b); }
  __forceinline const Color operator*(const float& a, const Color& b) { return Color(a*b.r,a*b.g,a*b.b); }
  __forceinline const Color operator*(const Color& a, const float& b) { return Color(a.r*b,a.g*b,a.b*b); }
  __forceinline const Color operator*(const Color& a, const Color& b) { return Color(a.r*b.r,a.g*b.g,a.b*b.b); }
  __forceinline const Color operator/(const Color& a, const Color& b) { return Color(a.r/b.r,a.g/b.g,a.b/b.b); }
  __forceinline const Color operator/(const Color& a, const float& b) { return Color(a.r/b,a.g/b,a.b/b); }

  __forceinline const Color min(const Color a, const Color b) { return Color(min(a.r,b.r),min(a.g,b.g),min(a.b,b.b)); }
  __forceinline const Color max(const Color a, const Color b) { return Color(max(a.r,b.r),max(a.g,b.g),max(a.b,b.b)); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Color operator+=(Color& a, const Color& b) { return a = a + b; }
  __forceinline const Color operator-=(Color& a, const Color& b) { return a = a - b; }
  __forceinline const Color operator*=(Color& a, const Color& b) { return a = a * b; }
  __forceinline const Color operator/=(Color& a, const Color& b) { return a = a / b; }
  __forceinline const Color operator*=(Color& a, const float& b) { return a = a * b; }
  __forceinline const Color operator/=(Color& a, const float& b) { return a = a / b; }

  ////////////////////////////////////////////////////////////////////////////////
  /// Reduction Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline float reduce_add(const Color& a) { return a.r+a.g+a.b; }
  __forceinline float reduce_mul(const Color& a) { return a.r*a.g*a.b; }
  __forceinline float reduce_min(const Color& a) { return min(a.r,a.g,a.b); }
  __forceinline float reduce_max(const Color& a) { return max(a.r,a.g,a.b); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline bool operator ==(const Color& a, const Color& b) { return a.r == b.r && a.g == b.g && a.b == b.b; }
  __forceinline bool operator !=(const Color& a, const Color& b) { return a.r != b.r || a.g != b.g || a.b != b.b; }
  __forceinline bool operator < (const Color& a, const Color& b ) {
    if (a.r != b.r) return a.r < b.r;
    if (a.g != b.g) return a.g < b.g;
    if (a.b != b.b) return a.b < b.b;
    return false;
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Select
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline Color select ( bool s, const Color& t, const Color& f ){ 
    return Color(select(s,t.r,f.r),select(s,t.g,f.g),select(s,t.b,f.b)); 
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Special Operators
  ////////////////////////////////////////////////////////////////////////////////

  /*! computes luminance of a color */
  __forceinline float luminance (const Color& a) { return 0.212671f*a.r + 0.715160f*a.g + 0.072169f*a.b; }

  __forceinline Color exp (const Color& a) { return Color(exp(a.r),exp(a.g),exp(a.b)); }
  __forceinline Color log (const Color& a) { return Color(log(a.r),log(a.g),log(a.b)); }
  __forceinline Color pow (const Color& a, float e) { return exp(log(max(Color(1E-10f),a))*e); }

  /*! output operator */
  inline std::ostream& operator<<(std::ostream& cout, const Color& a) {
    return cout << "(" << a.r << ", " << a.g << ", " << a.b << ")";
  }
}

#endif

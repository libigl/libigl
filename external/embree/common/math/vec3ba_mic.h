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

#include "math.h"

namespace embree
{
  ////////////////////////////////////////////////////////////////////////////////
  /// MIC Vec3ba Type
  ////////////////////////////////////////////////////////////////////////////////

  /*! 3-wide bool type emulated with 16-wide bit vectors. */
  struct Vec3ba 
  {
    enum { N = 3 };
    __mmask16 v; 

    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vec3ba( ) {}
    __forceinline Vec3ba           ( const __mmask16 a ) : v(a) {}
    __forceinline Vec3ba           ( const Vec3ba& a  ) : v(a.v) {}
    __forceinline Vec3ba& operator =(const Vec3ba& a) { v = a.v; return *this; }

    __forceinline explicit Vec3ba ( bool x                 ) : v(x ? 0xF : 0x0) {}
    __forceinline explicit Vec3ba ( bool x, bool y, bool z ) : v(1*__mmask16(x)+2*__mmask16(y)+4*__mmask16(z)) {}

    __forceinline operator const __mmask16&( void ) const { return v; }
    __forceinline operator       __mmask16&( void )       { return v; }

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vec3ba( FalseTy ) : v(0x0) {}
    __forceinline Vec3ba( TrueTy  ) : v(0xF) {}

    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline bool operator []( const size_t index ) const { assert(index < 3); return (v >> index) & 1; }
    //__forceinline       int& operator []( const size_t index )       { assert(index < 3); return (&x)[index]; }
  };


  ////////////////////////////////////////////////////////////////////////////////
  /// Unary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vec3ba operator !( const Vec3ba& a ) { return _mm512_knot(a); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vec3ba operator &( const Vec3ba& a, const Vec3ba& b ) { return _mm512_kand(a.v, b.v); }
  __forceinline const Vec3ba operator |( const Vec3ba& a, const Vec3ba& b ) { return _mm512_kor (a.v, b.v); }
  __forceinline const Vec3ba operator ^( const Vec3ba& a, const Vec3ba& b ) { return _mm512_kxor(a.v, b.v); }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline const Vec3ba operator &=( Vec3ba& a, const Vec3ba& b ) { return a = a & b; }
  __forceinline const Vec3ba operator |=( Vec3ba& a, const Vec3ba& b ) { return a = a | b; }
  __forceinline const Vec3ba operator ^=( Vec3ba& a, const Vec3ba& b ) { return a = a ^ b; }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators + Select
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline bool operator ==( const Vec3ba& a, const Vec3ba& b ) { return _mm512_kxnor(a, b); }
  __forceinline bool operator !=( const Vec3ba& a, const Vec3ba& b ) { return _mm512_kxor (a, b); }
  __forceinline bool operator < ( const Vec3ba& a, const Vec3ba& b ) {
    if (a[0] != b[0]) return a[0] < b[0];
    if (a[1] != b[1]) return a[1] < b[1];
    if (a[2] != b[2]) return a[2] < b[2];
    return false;
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Reduction Operations
  ////////////////////////////////////////////////////////////////////////////////
    
  __forceinline bool reduce_and( const Vec3ba& a ) { return (a.v & 0x7) == 0x7; }
  __forceinline bool reduce_or ( const Vec3ba& a ) { return (a.v & 0x7) != 0x0; }
  
  __forceinline int all(const Vec3ba& a)  { return (a.v & 0x7) == 0x7; }
  __forceinline int any(const Vec3ba& a)  { return (a.v & 0x7) != 0x0; }
  __forceinline int none(const Vec3ba& a) { return (a.v & 0x7) == 0x0; }

  ////////////////////////////////////////////////////////////////////////////////
  /// Output Operators
  ////////////////////////////////////////////////////////////////////////////////

  inline std::ostream& operator<<(std::ostream& cout, const Vec3ba& a) {
    return cout << "(" << (a[0] ? "1" : "0") << ", " << (a[1] ? "1" : "0") << ", " << (a[2] ? "1" : "0") << ")";
  }
}

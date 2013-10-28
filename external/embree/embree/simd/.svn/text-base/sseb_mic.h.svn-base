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

#ifndef __EMBREE_SSEB_MIC_H__
#define __EMBREE_SSEB_MIC_H__

namespace embree
{
  /* memory representation as bit vector */
  struct sseb_m 
  {
    enum { size = 4 };  // number of SIMD elements
    __mmask16 v; 
    
    __forceinline sseb_m ( ) {}
    __forceinline sseb_m ( bool  a                            ) : v(a ? 0xF : 0x0) {}
    __forceinline sseb_m ( bool  a, bool  b                   ) : v(1*__mmask16(a)+2*__mmask16(a)+4*__mmask16(b)+8*__mmask16(b)) {}
    __forceinline sseb_m ( bool  a, bool  b, bool  c, bool  d ) : v(1*__mmask16(a)+2*__mmask16(b)+4*__mmask16(c)+8*__mmask16(d)) {}

    __forceinline sseb_m           ( const __mmask16& other ) { v = other; }
    __forceinline sseb_m& operator=( const __mmask16& other ) { v = other; return *this; }
  };
  
  /*! 4-wide bool type emulated with 16-wide bit vectors. */
  struct sseb_t 
  {
    enum { size = 4 };  // number of SIMD elements
    __mmask16 v; 
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline sseb_t           ( ) {}
    __forceinline sseb_t           ( const sseb_t& other ) { v = other.v; }
    __forceinline sseb_t& operator=( const sseb_t& other ) { v = other.v; return *this; }
      
    __forceinline sseb_t           ( const sseb_m& other ) { v = other.v; }
    __forceinline sseb_t& operator=( const sseb_m& other ) {  v = other.v; return *this; }
      
    __forceinline sseb_t( const __mmask16 a ) : v(a) {}
    __forceinline operator const __mmask16&( void ) const { return v; }
    
    __forceinline sseb_t           ( bool  a                            ) : v(a ? 0xF : 0x0) {}
    __forceinline sseb_t           ( bool  a, bool  b                   ) : v(1*__mmask16(a)+2*__mmask16(a)+4*__mmask16(b)+8*__mmask16(b)) {}
    __forceinline sseb_t           ( bool  a, bool  b, bool  c, bool  d ) : v(1*__mmask16(a)+2*__mmask16(b)+4*__mmask16(c)+8*__mmask16(d)) {}
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline sseb_t( FalseTy ) : v(0x0) {}
    __forceinline sseb_t( TrueTy  ) : v(0xF) {}
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline bool   operator []( const size_t i ) const { assert(i < 4); return (v >> i) & 1; }
    //__forceinline int32& operator []( const size_t i )       { assert(i < 4); return v[i]; }
  };
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Unary Operators
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline const sseb_t operator !( const sseb_t& a ) { return _mm512_knot(a); }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline const sseb_t operator &( const sseb_t& a, const sseb_t& b ) { return _mm512_kand(a, b); }
  __forceinline const sseb_t operator |( const sseb_t& a, const sseb_t& b ) { return _mm512_kor (a, b); }
  __forceinline const sseb_t operator ^( const sseb_t& a, const sseb_t& b ) { return _mm512_kxor(a, b); }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline const sseb_t operator &=( sseb_t& a, const sseb_t& b ) { return a = a & b; }
  __forceinline const sseb_t operator |=( sseb_t& a, const sseb_t& b ) { return a = a | b; }
  __forceinline const sseb_t operator ^=( sseb_t& a, const sseb_t& b ) { return a = a ^ b; }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators + Select
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline const sseb_t operator !=( const sseb_t& a, const sseb_t& b ) { return _mm512_kxor (a, b); }
  __forceinline const sseb_t operator ==( const sseb_t& a, const sseb_t& b ) { return _mm512_kxnor(a, b); }
  __forceinline const sseb_t select( const sseb_t& m, const sseb_t& a, const sseb_t& b ) { return _mm512_kor(_mm512_kand(m, a), _mm512_kandn(m, b)); }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Reduction Operations
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline size_t popcnt( const sseb_t& a ) { return __popcnt(a); }
  
  __forceinline bool reduce_and( const sseb_t& a ) { return (a.v & 0xf) == 0xf; }
  __forceinline bool reduce_or ( const sseb_t& a ) { return (a.v & 0xf) != 0x0; }
  __forceinline bool all       ( const sseb_t& a ) { return (a.v & 0xf) == 0xf; }
  __forceinline bool any       ( const sseb_t& a ) { return (a.v & 0xf) != 0x0; }
  __forceinline bool none      ( const sseb_t& a ) { return (a.v & 0xf) == 0x0; }
  
  __forceinline size_t movemask( const sseb_t& a ) { return a.v & 0xf; }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Output Operators
  ////////////////////////////////////////////////////////////////////////////////
  
  inline std::ostream& operator<<(std::ostream& cout, const sseb_t& a) {
    return cout << "<" << a[0] << ", " << a[1] << ", " << a[2] << ", " << a[3] << ">";
  }
}

#endif

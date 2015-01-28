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
  struct __aligned(16) sseb_m 
  {
    enum { size = 4 };  // number of SIMD elements
    int32 v[4];
    
    __forceinline sseb_m ( ) {}

    __forceinline sseb_m ( bool  a ) { 
      v[0] = a ? -1 : 0;
      v[1] = a ? -1 : 0;
      v[2] = a ? -1 : 0;
      v[3] = a ? -1 : 0;
    }
    __forceinline sseb_m ( bool  a, bool  b ) {
      v[0] = a ? -1 : 0;
      v[1] = b ? -1 : 0;
      v[2] = a ? -1 : 0;
      v[3] = b ? -1 : 0;
    }
    __forceinline sseb_m ( bool  a, bool  b, bool  c, bool  d ) {
      v[0] = a ? -1 : 0;
      v[1] = b ? -1 : 0;
      v[2] = c ? -1 : 0;
      v[3] = d ? -1 : 0;
    }

    __forceinline sseb_m           ( const sseb_t& other );
    __forceinline sseb_m& operator=( const sseb_t& other );

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline sseb_m( FalseTy ) { v[0] = v[1] = v[2] = v[3] = 0; } 
    __forceinline sseb_m( TrueTy  ) { v[0] = v[1] = v[2] = v[3] = -1; } 
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline bool   operator []( const size_t i ) const { assert(i < 4); return v[i] != 0; }
    __forceinline int32& operator []( const size_t i )       { assert(i < 4); return v[i]; }
  };
  
  /*! 4-wide bool type emulated with 16-wide bit vectors. */
  struct sseb_t 
  {
    __mmask16 v; 
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline sseb_t           ( const __mmask16 a ) : v(a) {}
    __forceinline sseb_t           ( const sseb_t& other ) { v = other.v; }
    __forceinline sseb_t& operator=( const sseb_t& other ) { v = other.v; return *this; }
    __forceinline sseb_t           ( const sseb_m& other ) { 
      __m512i i = _mm512_extload_epi32(&other,_MM_UPCONV_EPI32_NONE,_MM_BROADCAST_4X16,_MM_HINT_NONE); 
      v = _mm512_cmp_epi32_mask(i,_mm512_set1_epi32(0),_MM_CMPINT_NE);
    }
  };
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline sseb_m::sseb_m ( const sseb_t& other ) { 
   __m512i i = _mm512_mask_blend_epi32(other.v, _mm512_set1_epi32(0), _mm512_set1_epi32(-1));
   _mm512_mask_extpackstorelo_epi32(this,0xf,i,_MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE); 
  }
  __forceinline sseb_m& sseb_m::operator=( const sseb_t& other ) { 
    __m512i i = _mm512_mask_blend_epi32(other.v, _mm512_set1_epi32(0), _mm512_set1_epi32(-1));
    _mm512_mask_extpackstorelo_epi32(this,0xf,i,_MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE); 
    return *this; 
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Unary Operators
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline const sseb_t operator !( const sseb_t& a ) { return _mm512_knot(a.v); }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline const sseb_t operator &( const sseb_t& a, const sseb_t& b ) { return _mm512_kand(a.v, b.v); }
  __forceinline const sseb_t operator |( const sseb_t& a, const sseb_t& b ) { return _mm512_kor (a.v, b.v); }
  __forceinline const sseb_t operator ^( const sseb_t& a, const sseb_t& b ) { return _mm512_kxor(a.v, b.v); }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline const sseb_t operator &=( sseb_m& a, const sseb_t& b ) { return a = a & b; }
  __forceinline const sseb_t operator |=( sseb_m& a, const sseb_t& b ) { return a = a | b; }
  __forceinline const sseb_t operator ^=( sseb_m& a, const sseb_t& b ) { return a = a ^ b; }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators + Select
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline const sseb_t operator !=( const sseb_t& a, const sseb_t& b ) { return _mm512_kxor (a.v, b.v); }
  __forceinline const sseb_t operator ==( const sseb_t& a, const sseb_t& b ) { return _mm512_kxnor(a.v, b.v); }
  __forceinline const sseb_t select( const sseb_t& m, const sseb_t& a, const sseb_t& b ) { return _mm512_kor(_mm512_kand(m.v, a.v), _mm512_kandn(m.v, b.v)); }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Reduction Operations
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline bool reduce_and( const sseb_t& a ) { return (a.v & 0xf) == 0xf; }
  __forceinline bool reduce_or ( const sseb_t& a ) { return (a.v & 0xf) != 0x0; }
  __forceinline bool all       ( const sseb_t& a ) { return (a.v & 0xf) == 0xf; }
  __forceinline bool any       ( const sseb_t& a ) { return (a.v & 0xf) != 0x0; }
  __forceinline bool none      ( const sseb_t& a ) { return (a.v & 0xf) == 0x0; }
  
  __forceinline size_t movemask( const sseb_t& a ) { return a.v & 0xf; }
  __forceinline size_t popcnt( const sseb_t& a ) { return __popcnt(a.v); }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Output Operators
  ////////////////////////////////////////////////////////////////////////////////
  
  inline std::ostream& operator<<(std::ostream& cout, const sseb_m& a) {
    return cout << "<" << a[0] << ", " << a[1] << ", " << a[2] << ", " << a[3] << ">";
  }
}

#endif

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
  /// MIC Vec3ia Type
  ////////////////////////////////////////////////////////////////////////////////

  struct Vec3ia_t;
  struct Vec3fa_t;

  /* 3 aligned ints as memory representation */
  struct __aligned(16) Vec3ia 
  {
    typedef int Scalar;
    enum { N = 3 }; 
    struct { int x,y,z; int a; };
        
    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline Vec3ia ( ) {}
    __forceinline Vec3ia ( int x               ) : x(x), y(x), z(x), a(0) {}
    __forceinline Vec3ia ( int x, int y, int z ) : x(x), y(y), z(z), a(0) {}

    __forceinline Vec3ia           ( const Vec3fa_t& other ); 
    __forceinline Vec3ia           ( const Vec3ia_t& other ); 
    __forceinline Vec3ia& operator=( const Vec3ia_t& other );
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vec3ia( ZeroTy   ) : x(0),       y(0),       z(0),       a(0) {}
    __forceinline Vec3ia( OneTy    ) : x(1),       y(1),       z(1),       a(1) {}
    __forceinline Vec3ia( PosInfTy ) : x(pos_inf), y(pos_inf), z(pos_inf), a(pos_inf) {}
    __forceinline Vec3ia( NegInfTy ) : x(neg_inf), y(neg_inf), z(neg_inf), a(neg_inf) {}

    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline const int& operator []( const size_t index ) const { assert(index < 3); return (&x)[index]; }
    __forceinline       int& operator []( const size_t index )       { assert(index < 3); return (&x)[index]; }
  };
  
  /*! 3-wide vectors emulated with 16-wide vectors. */
  struct __aligned(64) Vec3ia_t 
  {
    __m512i m512; 
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vec3ia_t            ( const __m512i a ) : m512(a) {}
    __forceinline Vec3ia_t            ( const Vec3ia_t& other ) { m512 = other.m512; }
    __forceinline Vec3ia_t& operator =( const Vec3ia_t& other ) { m512 = other.m512; return *this; }
    
    __forceinline operator const __m512i&( void ) const { return m512; }
    __forceinline operator       __m512i&( void )       { return m512; }

  public:
    __forceinline explicit Vec3ia_t ( int a ) 
      : m512(_mm512_set1_epi32(a)) {}
    
    __forceinline Vec3ia_t( const Vec3ia& other ) { 
      m512 = _mm512_extload_epi32(&other,_MM_UPCONV_EPI32_NONE,_MM_BROADCAST_4X16,_MM_HINT_NONE); 
    }
  };

  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline Vec3ia::Vec3ia ( const Vec3ia_t& other ) { 
    _mm512_mask_extpackstorelo_epi32(this,0xf,other.m512,_MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE); 
  }
  
  __forceinline Vec3ia& Vec3ia::operator=( const Vec3ia_t& other ) { 
    _mm512_mask_extpackstorelo_epi32(this,0xf,other.m512,_MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE); 
    return *this; 
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Unary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vec3ia_t operator +( const Vec3ia_t& a ) { return a; }
  __forceinline const Vec3ia_t operator -( const Vec3ia_t& a ) { return _mm512_sub_epi32(_mm512_setzero_epi32(), a.m512); }
  //__forceinline const Vec3ia_t abs       ( const Vec3ia_t& a ) { return _mm512_abs_epi32(a.m512); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vec3ia_t operator +( const Vec3ia_t& a, const Vec3ia_t& b ) { return _mm512_add_epi32(a.m512, b.m512); }
  __forceinline const Vec3ia_t operator +( const Vec3ia_t& a, const int      b ) { return a+Vec3ia_t(b); }
  __forceinline const Vec3ia_t operator +( const int      a, const Vec3ia_t& b ) { return Vec3ia_t(a)+b; }

  __forceinline const Vec3ia_t operator -( const Vec3ia_t& a, const Vec3ia_t& b ) { return _mm512_sub_epi32(a.m512, b.m512); }
  __forceinline const Vec3ia_t operator -( const Vec3ia_t& a, const int      b ) { return a-Vec3ia_t(b); }
  __forceinline const Vec3ia_t operator -( const int      a, const Vec3ia_t& b ) { return Vec3ia_t(a)-b; }

  __forceinline const Vec3ia_t operator *( const Vec3ia_t& a, const Vec3ia_t& b ) { return _mm512_mullo_epi32(a.m512, b.m512); }
  __forceinline const Vec3ia_t operator *( const Vec3ia_t& a, const int      b ) { return a * Vec3ia_t(b); }
  __forceinline const Vec3ia_t operator *( const int      a, const Vec3ia_t& b ) { return Vec3ia_t(a) * b; }

  __forceinline const Vec3ia_t operator &( const Vec3ia_t& a, const Vec3ia_t& b ) { return _mm512_and_epi32(a.m512, b.m512); }
  __forceinline const Vec3ia_t operator &( const Vec3ia_t& a, const int      b ) { return a & Vec3ia_t(b); }
  __forceinline const Vec3ia_t operator &( const int      a, const Vec3ia_t& b ) { return Vec3ia_t(a) & b; }

  __forceinline const Vec3ia_t operator |( const Vec3ia_t& a, const Vec3ia_t& b ) { return _mm512_or_epi32(a.m512, b.m512); }
  __forceinline const Vec3ia_t operator |( const Vec3ia_t& a, const int      b ) { return a | Vec3ia_t(b); }
  __forceinline const Vec3ia_t operator |( const int      a, const Vec3ia_t& b ) { return Vec3ia_t(a) | b; }

  __forceinline const Vec3ia_t operator ^( const Vec3ia_t& a, const Vec3ia_t& b ) { return _mm512_xor_epi32(a.m512, b.m512); }
  __forceinline const Vec3ia_t operator ^( const Vec3ia_t& a, const int      b ) { return a ^ Vec3ia_t(b); }
  __forceinline const Vec3ia_t operator ^( const int      a, const Vec3ia_t& b ) { return Vec3ia_t(a) ^ b; }

  __forceinline const Vec3ia_t operator <<( const Vec3ia_t& a, const int n ) { return _mm512_slli_epi32(a.m512, n); }
  __forceinline const Vec3ia_t operator >>( const Vec3ia_t& a, const int n ) { return _mm512_srai_epi32(a.m512, n); }

  __forceinline const Vec3ia_t sra ( const Vec3ia_t& a, const int b ) { return _mm512_srai_epi32(a.m512, b); }
  __forceinline const Vec3ia_t srl ( const Vec3ia_t& a, const int b ) { return _mm512_srli_epi32(a.m512, b); }

  __forceinline const Vec3ia_t min( const Vec3ia_t& a, const Vec3ia_t& b ) { return _mm512_min_epi32(a.m512,b.m512); }
  __forceinline const Vec3ia_t max( const Vec3ia_t& a, const Vec3ia_t& b ) { return _mm512_max_epi32(a.m512,b.m512); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline Vec3ia& operator +=( Vec3ia& a, const Vec3ia_t& b ) { return a = a + b; }
  __forceinline Vec3ia& operator +=( Vec3ia& a, const int32&   b ) { return a = a + b; }
  
  __forceinline Vec3ia& operator -=( Vec3ia& a, const Vec3ia_t& b ) { return a = a - b; }
  __forceinline Vec3ia& operator -=( Vec3ia& a, const int32&   b ) { return a = a - b; }
  
  __forceinline Vec3ia& operator *=( Vec3ia& a, const Vec3ia_t& b ) { return a = a * b; }
  __forceinline Vec3ia& operator *=( Vec3ia& a, const int32&   b ) { return a = a * b; }
  
  __forceinline Vec3ia& operator &=( Vec3ia& a, const Vec3ia_t& b ) { return a = a & b; }
  __forceinline Vec3ia& operator &=( Vec3ia& a, const int32&   b ) { return a = a & b; }
  
  __forceinline Vec3ia& operator |=( Vec3ia& a, const Vec3ia_t& b ) { return a = a | b; }
  __forceinline Vec3ia& operator |=( Vec3ia& a, const int32&   b ) { return a = a | b; }
  
  __forceinline Vec3ia& operator <<=( Vec3ia& a, const int32&  b ) { return a = a << b; }
  __forceinline Vec3ia& operator >>=( Vec3ia& a, const int32&  b ) { return a = a >> b; }

  ////////////////////////////////////////////////////////////////////////////////
  /// Reductions
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline int reduce_add(const Vec3ia_t& a) { 
    __m512i b = _mm512_mask_blend_epi32(0x8,a.m512,_mm512_set1_epi32(0));
    __m512i c = _mm512_add_epi32(b, _mm512_swizzle_epi32(b, _MM_SWIZ_REG_BADC));
    __m512i d = _mm512_add_epi32(c, _mm512_swizzle_epi32(c, _MM_SWIZ_REG_CDAB));
    return *(int*)&d;
  } 

  __forceinline int reduce_mul(const Vec3ia_t& a) { 
    __m512i b = _mm512_mask_blend_epi32(0x8,a.m512,_mm512_set1_epi32(1));
    __m512i c = _mm512_mullo_epi32(b, _mm512_swizzle_epi32(b, _MM_SWIZ_REG_BADC));
    __m512i d = _mm512_mullo_epi32(c, _mm512_swizzle_epi32(c, _MM_SWIZ_REG_CDAB));
    return *(int*)&d;
  }

  __forceinline int reduce_min(const Vec3ia_t& a) { 
    __m512i b = _mm512_mask_blend_epi32(0x8,a.m512,_mm512_set1_epi32(pos_inf));
    __m512i c = _mm512_min_epi32(b, _mm512_swizzle_epi32(b, _MM_SWIZ_REG_BADC));
    __m512i d = _mm512_min_epi32(c, _mm512_swizzle_epi32(c, _MM_SWIZ_REG_CDAB));
    return *(int*)&d;
  } 

  __forceinline int reduce_max(const Vec3ia_t& a) { 
    __m512i b = _mm512_mask_blend_epi32(0x8,a.m512,_mm512_set1_epi32(neg_inf));
    __m512i c = _mm512_max_epi32(b, _mm512_swizzle_epi32(b, _MM_SWIZ_REG_BADC));
    __m512i d = _mm512_max_epi32(c, _mm512_swizzle_epi32(c, _MM_SWIZ_REG_CDAB));
    return *(int*)&d;
  } 

  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline bool operator ==( const Vec3ia_t& a, const Vec3ia_t& b ) { return (_mm512_cmp_epi32_mask(a,b,_MM_CMPINT_EQ) & 7) == 7; }
  __forceinline bool operator !=( const Vec3ia_t& a, const Vec3ia_t& b ) { return (_mm512_cmp_epi32_mask(a,b,_MM_CMPINT_NE) & 7) != 0; }

  ////////////////////////////////////////////////////////////////////////////////
  /// Select
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vec3ia_t select( bool s, const Vec3ia_t& t, const Vec3ia_t& f ) {
    return _mm512_mask_blend_epi32(s ? 0xF : 0x0, f, t);
  }

  __forceinline const Vec3ia_t select( const Vec3ba& s, const Vec3ia_t& t, const Vec3ia_t& f ) {
    return _mm512_mask_blend_epi32(s, f, t);
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Output Operators
  ////////////////////////////////////////////////////////////////////////////////

  inline std::ostream& operator<<(std::ostream& cout, const Vec3ia& a) {
    return cout << "(" << a.x << ", " << a.y << ", " << a.z << ")";
  }
}

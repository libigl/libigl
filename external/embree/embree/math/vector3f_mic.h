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

#ifndef __EMBREE_VEC3F_MIC_H__
#define __EMBREE_VEC3F_MIC_H__

#include "math.h"

namespace embree
{
  ////////////////////////////////////////////////////////////////////////////////
  /// MIC Vector3f Type
  ////////////////////////////////////////////////////////////////////////////////

  struct Vector3f;

  /* 3 aligned floats as memory representation */
  struct __align(16) Vec3fa 
  {
    typedef float Scalar;
    enum { N = 3 };
    float x,y,z; int a;
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline Vec3fa ( ) {}
    __forceinline Vec3fa ( float x                   ) : x(x), y(x), z(x), a(0) {}
    __forceinline Vec3fa ( float x, float y, float z ) : x(x), y(y), z(z), a(0) {}

    __forceinline Vec3fa ( const Vector3f& other ); 
    __forceinline Vec3fa& operator=( const Vector3f& other );

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    /* __forceinline Vec3fa( ZeroTy   ) : x(0.0f),    y(0.0f),    z(0.0f),    a(0) {}
    __forceinline Vec3fa( OneTy    ) : x(1.0f),    y(1.0f),    z(1.0f),    a(0) {}
    __forceinline Vec3fa( PosInfTy ) : x(pos_inf), y(pos_inf), z(pos_inf), a(0) {}
    __forceinline Vec3fa( NegInfTy ) : x(neg_inf), y(neg_inf), z(neg_inf), a(0) {}*/

    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline const float& operator []( const size_t index ) const { assert(index < 3); return (&x)[index]; }
    __forceinline       float& operator []( const size_t index )       { assert(index < 3); return (&x)[index]; }
  };
  
  /*! 3-wide vectors emulated with 16-wide vectors. */
  struct __align(64) Vector3f 
  {
    typedef float Scalar;
    enum { N = 3 };
    union {
      __m512 m512; 
      struct { float x,y,z; int a; };
    };

    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vector3f( ) {}
    __forceinline Vector3f( const __m512 a ) : m512(a) {}

    __forceinline Vector3f            ( const Vector3f& other ) { m512 = other.m512; }
    __forceinline Vector3f& operator =( const Vector3f& other ) { m512 = other.m512; return *this; }

    __forceinline Vector3f            ( const Vec3fa& other ) { m512 = _mm512_extload_ps(&other,_MM_UPCONV_PS_NONE,_MM_BROADCAST_4X16,_MM_HINT_NONE); }
    __forceinline Vector3f& operator =( const Vec3fa& other ) { m512 = _mm512_extload_ps(&other,_MM_UPCONV_PS_NONE,_MM_BROADCAST_4X16,_MM_HINT_NONE); return *this; }

    __forceinline explicit Vector3f( const float a ) : m512(_mm512_set1_ps(a)) {}
    __forceinline explicit Vector3f( const float x, const float y, const float z) : m512(Vector3f(Vec3fa(x,y,z))) {}
    __forceinline explicit Vector3f( const float* const a, const ssize_t stride = 1 ) : m512(Vector3f(Vec3fa(a[0*stride],a[1*stride],a[2*stride]))) {}

    __forceinline explicit Vector3f( const __m512i a ) : m512(_mm512_cvtfxpnt_round_adjustepi32_ps(a,_MM_ROUND_MODE_NEAREST,_MM_EXPADJ_NONE)) {}

    __forceinline operator const __m512&( void ) const { return m512; }
    __forceinline operator       __m512&( void )       { return m512; }

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vector3f( ZeroTy   ) : m512(_mm512_setzero_ps()) {}
    __forceinline Vector3f( OneTy    ) : m512(_mm512_set1_ps(1.0f)) {}
    __forceinline Vector3f( PosInfTy ) : m512(_mm512_set1_ps(pos_inf)) {}
    __forceinline Vector3f( NegInfTy ) : m512(_mm512_set1_ps(neg_inf)) {}

    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline const float& operator []( const size_t index ) const { assert(index < 3); return (&x)[index]; }
    __forceinline       float& operator []( const size_t index )       { assert(index < 3); return (&x)[index]; }
  };

  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline Vec3fa::Vec3fa ( const Vector3f& other ) { 
    assert((size_t)this % 16 == 0); 
    _mm512_mask_extpackstorelo_ps(this,0xf,other.m512,_MM_DOWNCONV_PS_NONE, _MM_HINT_NONE); 
  }
  
  __forceinline Vec3fa& Vec3fa::operator=( const Vector3f& other ) { 
    assert((size_t)this % 16 == 0); 
    _mm512_mask_extpackstorelo_ps(this,0xf,other.m512,_MM_DOWNCONV_PS_NONE, _MM_HINT_NONE); 
    return *this; 
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Unary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vector3f operator +( const Vector3f& a ) { return a; }
  __forceinline const Vector3f operator -( const Vector3f& a ) { return _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(a.m512), _mm512_set1_epi32(0x80000000))); }
  __forceinline const Vector3f abs       ( const Vector3f& a ) { return _mm512_gmaxabs_ps(a.m512,a.m512); }
  __forceinline const Vector3f sign      ( const Vector3f& a ) { return _mm512_mask_blend_ps(_mm512_cmp_ps_mask(a,_mm512_setzero_ps(),_MM_CMPINT_LT), _mm512_set1_ps(1.0f), _mm512_set1_ps(-1.0f)); }

  __forceinline const Vector3f rcp  ( const Vector3f& a ) { return  _mm512_rcp23_ps(a); }
  __forceinline const Vector3f sqr  ( const Vector3f& a ) { return _mm512_mul_ps(a,a); }
  __forceinline const Vector3f sqrt ( const Vector3f& a ) { return _mm512_sqrt_ps(a.m512); }
  __forceinline const Vector3f rsqrt( const Vector3f& a ) { return _mm512_invsqrt_ps(a.m512); }
  __forceinline const Vector3f zero_fix( const Vector3f& a ) { return _mm512_mask_blend_ps(_mm512_cmp_ps_mask(a.m512,_mm512_setzero_ps(),_MM_CMPINT_EQ), a, _mm512_set1_ps(1E-10f)); }
  __forceinline const Vector3f rcp_safe(const Vector3f& a) { return rcp(zero_fix(a)); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vector3f operator +( const Vector3f& a, const Vector3f& b ) { return _mm512_add_ps(a.m512, b.m512); }
  __forceinline const Vector3f operator -( const Vector3f& a, const Vector3f& b ) { return _mm512_sub_ps(a.m512, b.m512); }
  __forceinline const Vector3f operator *( const Vector3f& a, const Vector3f& b ) { return _mm512_mul_ps(a.m512, b.m512); }
  __forceinline const Vector3f operator *( const Vector3f& a, const float    b ) { return a * Vector3f(b); }
  __forceinline const Vector3f operator *( const float    a, const Vector3f& b ) { return Vector3f(a) * b; }
  __forceinline const Vector3f operator /( const Vector3f& a, const Vector3f& b ) { return _mm512_div_ps(a.m512,b.m512); }
  __forceinline const Vector3f operator /( const Vector3f& a, const float    b ) { return _mm512_div_ps(a.m512,_mm512_set1_ps(b)); }
  __forceinline const Vector3f operator /( const    float a, const Vector3f& b ) { return _mm512_div_ps(_mm512_set1_ps(a),b.m512); }

  __forceinline const Vector3f min( const Vector3f& a, const Vector3f& b ) { return _mm512_gmin_ps(a.m512,b.m512); }
  __forceinline const Vector3f max( const Vector3f& a, const Vector3f& b ) { return _mm512_gmax_ps(a.m512,b.m512); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline Vector3f& operator +=( Vector3f& a, const Vector3f& b ) { return a = a + b; }
  __forceinline Vector3f& operator -=( Vector3f& a, const Vector3f& b ) { return a = a - b; }
  __forceinline Vector3f& operator *=( Vector3f& a, const float    b ) { return a = a * b; }
  __forceinline Vector3f& operator /=( Vector3f& a, const float    b ) { return a = a / b; }

  ////////////////////////////////////////////////////////////////////////////////
  /// Reductions
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline float reduce_add(const Vector3f& a) { 
    __m512 b = _mm512_mask_blend_ps(0x8,a.m512,_mm512_set1_ps(0.0f));
    __m512 c = _mm512_add_ps(b, _mm512_swizzle_ps(b, _MM_SWIZ_REG_BADC));
    __m512 d = _mm512_add_ps(c, _mm512_swizzle_ps(c, _MM_SWIZ_REG_CDAB));
    return _mm512_cvtss_f32(d);
  } 

  __forceinline float reduce_mul(const Vector3f& a) { 
    __m512 b = _mm512_mask_blend_ps(0x8,a.m512,_mm512_set1_ps(1.0f));
    __m512 c = _mm512_mul_ps(b, _mm512_swizzle_ps(b, _MM_SWIZ_REG_BADC));
    __m512 d = _mm512_mul_ps(c, _mm512_swizzle_ps(c, _MM_SWIZ_REG_CDAB));
    return _mm512_cvtss_f32(d);
  } 

  __forceinline float reduce_min(const Vector3f& a) { 
    __m512 b = _mm512_mask_blend_ps(0x8,a.m512,_mm512_set1_ps(pos_inf));
    __m512 c = _mm512_gmin_ps(b, _mm512_swizzle_ps(b, _MM_SWIZ_REG_BADC));
    __m512 d = _mm512_gmin_ps(c, _mm512_swizzle_ps(c, _MM_SWIZ_REG_CDAB));
    return _mm512_cvtss_f32(d);
  } 

  __forceinline float reduce_max(const Vector3f& a) { 
    __m512 b = _mm512_mask_blend_ps(0x8,a.m512,_mm512_set1_ps(neg_inf));
    __m512 c = _mm512_gmax_ps(b, _mm512_swizzle_ps(b, _MM_SWIZ_REG_BADC));
    __m512 d = _mm512_gmax_ps(c, _mm512_swizzle_ps(c, _MM_SWIZ_REG_CDAB));
    return _mm512_cvtss_f32(d);
  } 

  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline bool operator ==( const Vector3f& a, const Vector3f& b ) { return (_mm512_cmp_ps_mask(a,b,_MM_CMPINT_EQ) & 7) == 7; }
  __forceinline bool operator !=( const Vector3f& a, const Vector3f& b ) { return (_mm512_cmp_ps_mask(a,b,_MM_CMPINT_NE) & 7) != 0; }
  __forceinline bool operator < ( const Vector3f& a, const Vector3f& b ) {
    if (a.x != b.x) return a.x < b.x;
    if (a.y != b.y) return a.y < b.y;
    if (a.z != b.z) return a.z < b.z;
    return false;
  }

  __forceinline Vector3b eq_mask( const Vector3f& a, const Vector3f& b ) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_EQ); }
  __forceinline Vector3b neq_mask(const Vector3f& a, const Vector3f& b ) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_NE); }
  __forceinline Vector3b lt_mask( const Vector3f& a, const Vector3f& b ) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_LT); }
  __forceinline Vector3b le_mask( const Vector3f& a, const Vector3f& b ) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_LE); }
  __forceinline Vector3b gt_mask( const Vector3f& a, const Vector3f& b ) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_GT); }
  __forceinline Vector3b ge_mask( const Vector3f& a, const Vector3f& b ) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_GE); }

 ////////////////////////////////////////////////////////////////////////////////
  /// Euclidian Space Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline float dot ( const Vector3f& a, const Vector3f& b ) {
    return reduce_add(a*b);
  }

  __forceinline Vector3f cross ( const Vector3f& a, const Vector3f& b ) 
  {
    __m512 c = _mm512_mul_ps(b,_mm512_swizzle_ps(a,_MM_SWIZ_REG_DACB));
    c = _mm512_fmsub_ps(_mm512_swizzle_ps(b,_MM_SWIZ_REG_DACB),a,c);
    c = _mm512_swizzle_ps(c,_MM_SWIZ_REG_DACB);
    return c;
  }

  __forceinline float   length   ( const Vector3f& a )                   { return sqrt(dot(a,a)); }
  __forceinline Vector3f normalize( const Vector3f& a )                   { return a*rsqrt(dot(a,a)); }
  __forceinline float   distance ( const Vector3f& a, const Vector3f& b ) { return length(a-b); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Select
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vector3f select( bool s, const Vector3f& t, const Vector3f& f ) {
    return _mm512_mask_blend_ps(s ? 0xF : 0x0, f, t);
  }

  __forceinline const Vector3f select( const Vector3b& s, const Vector3f& t, const Vector3f& f ) {
    return _mm512_mask_blend_ps(s, f, t);
  }

  __forceinline int maxDim ( const Vector3f& a ) 
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

  inline std::ostream& operator<<(std::ostream& cout, const Vector3f& a) {
    return cout << "(" << a.x << ", " << a.y << ", " << a.z << ")";
  }
}

#endif

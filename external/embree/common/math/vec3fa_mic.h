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
#include "simd/mic.h"

namespace embree
{
  ////////////////////////////////////////////////////////////////////////////////
  /// MIC Vec3fa Type
  ////////////////////////////////////////////////////////////////////////////////

  struct Vec3fa_t;

  /* 3 aligned floats as memory representation */
  struct __aligned(16) Vec3fa 
  {
    typedef float Scalar;
    enum { N = 3 };
    float x,y,z; union { int a; unsigned u; float w; };
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline Vec3fa () {}
    __forceinline Vec3fa ( float x                   ) : x(x), y(x), z(x), a(0) {}
    __forceinline Vec3fa ( float x, float y, float z ) : x(x), y(y), z(z), a(0) {}
    __forceinline Vec3fa ( const Vec3f & o           ) : x(o.x), y(o.y), z(o.z), a(0) {}
    __forceinline Vec3fa ( const Vec3fa& other, const int other_a ) : x(other.x), y(other.y), z(other.z), a(other_a) {}
    
    __forceinline Vec3fa           ( const Vec3ia_t& other );
    __forceinline Vec3fa           ( const Vec3fa_t& other ); 
    __forceinline Vec3fa& operator=( const Vec3fa_t& other );

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline Vec3fa( ZeroTy   ) : x(0.0f),    y(0.0f),    z(0.0f),    a(0) {}
    __forceinline Vec3fa( OneTy    ) : x(1.0f),    y(1.0f),    z(1.0f),    a(0) {}
    __forceinline Vec3fa( PosInfTy ) : x(pos_inf), y(pos_inf), z(pos_inf), a(0) {}
    __forceinline Vec3fa( NegInfTy ) : x(neg_inf), y(neg_inf), z(neg_inf), a(0) {}

    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline const float& operator []( const size_t index ) const { assert(index < 3); return (&x)[index]; }
    __forceinline       float& operator []( const size_t index )       { assert(index < 3); return (&x)[index]; }
  };
  
  /*! 3-wide vectors emulated with 16-wide vectors. */
  struct __aligned(64) Vec3fa_t 
  {
    __m512 m512; 
        
    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////
    __forceinline Vec3fa_t( ) { m512 = _mm512_undefined(); }

    __forceinline Vec3fa_t( const __m512 a ) : m512(a) {}
    __forceinline Vec3fa_t            ( const Vec3fa_t& other ) { m512 = other.m512; }
    __forceinline Vec3fa_t& operator =( const Vec3fa_t& other ) { m512 = other.m512; return *this; }
    
    __forceinline operator const __m512&( void ) const { return m512; }
    __forceinline operator       __m512&( void )       { return m512; }

    __forceinline operator const mic_f( void ) const { return mic_f(m512); }
    __forceinline operator       mic_f( void )       { return mic_f(m512); }

  public:
    __forceinline explicit Vec3fa_t ( float a ) 
      : m512(_mm512_set1_ps(a)) {}
    
    __forceinline Vec3fa_t( const Vec3fa& other ) { 
      m512 = _mm512_extload_ps(&other,_MM_UPCONV_PS_NONE,_MM_BROADCAST_4X16,_MM_HINT_NONE); 
    }
  };

  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline Vec3ia::Vec3ia ( const Vec3fa_t& other ) {
    __m512i i = _mm512_cvtfxpnt_round_adjustps_epi32(other,_MM_ROUND_MODE_NEAREST,_MM_EXPADJ_NONE);
    _mm512_mask_extpackstorelo_epi32(this,0xf,i,_MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE); 
  }

  __forceinline Vec3fa::Vec3fa( const Vec3ia_t& v ) {
    __m512 f = _mm512_cvtfxpnt_round_adjustepi32_ps(v,_MM_ROUND_MODE_NEAREST,_MM_EXPADJ_NONE);
    _mm512_mask_extpackstorelo_ps(this,0xf,f,_MM_DOWNCONV_PS_NONE, _MM_HINT_NONE); 
  }

  __forceinline Vec3fa::Vec3fa ( const Vec3fa_t& other ) { 
    _mm512_mask_extpackstorelo_ps(this,0xf,other.m512,_MM_DOWNCONV_PS_NONE, _MM_HINT_NONE); 
  }
  
  __forceinline Vec3fa& Vec3fa::operator=( const Vec3fa_t& other ) { 
    _mm512_mask_extpackstorelo_ps(this,0xf,other.m512,_MM_DOWNCONV_PS_NONE, _MM_HINT_NONE); 
    return *this; 
  }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Unary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vec3fa_t operator +( const Vec3fa_t& a ) { return a; }
  __forceinline const Vec3fa_t operator -( const Vec3fa_t& a ) { return _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(a.m512), _mm512_set1_epi32(0x80000000))); }
  __forceinline const Vec3fa_t abs       ( const Vec3fa_t& a ) { return _mm512_gmaxabs_ps(a.m512,a.m512); }
  __forceinline const Vec3fa_t sign      ( const Vec3fa_t& a ) { return _mm512_mask_blend_ps(_mm512_cmp_ps_mask(a,_mm512_setzero_ps(),_MM_CMPINT_LT), _mm512_set1_ps(1.0f), _mm512_set1_ps(-1.0f)); }

  __forceinline const Vec3fa_t rcp  ( const Vec3fa_t& a ) { return  _mm512_rcp23_ps(a); }
  __forceinline const Vec3fa_t sqr  ( const Vec3fa_t& a ) { return _mm512_mul_ps(a,a); }
  __forceinline const Vec3fa_t sqrt ( const Vec3fa_t& a ) { return _mm512_sqrt_ps(a.m512); }
  __forceinline const Vec3fa_t rsqrt( const Vec3fa_t& a ) { return _mm512_invsqrt_ps(a.m512); }
  __forceinline const Vec3fa_t zero_fix( const Vec3fa_t& a ) { return _mm512_mask_blend_ps(_mm512_cmp_ps_mask(a.m512,_mm512_setzero_ps(),_MM_CMPINT_EQ), a, _mm512_set1_ps(1E-10f)); }
  __forceinline const Vec3fa_t rcp_safe(const Vec3fa_t& a) { return rcp(zero_fix(a)); }

  __forceinline Vec3fa log ( const Vec3fa& a ) { 
    return Vec3fa(logf(a.x),logf(a.y),logf(a.z));
  }

  __forceinline Vec3fa exp ( const Vec3fa& a ) { 
    return Vec3fa(expf(a.x),expf(a.y),expf(a.z));
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vec3fa_t operator +( const Vec3fa_t& a, const Vec3fa_t& b ) { return _mm512_add_ps(a.m512, b.m512); }
  __forceinline const Vec3fa_t operator -( const Vec3fa_t& a, const Vec3fa_t& b ) { return _mm512_sub_ps(a.m512, b.m512); }
  __forceinline const Vec3fa_t operator *( const Vec3fa_t& a, const Vec3fa_t& b ) { return _mm512_mul_ps(a.m512, b.m512); }
  __forceinline const Vec3fa_t operator *( const Vec3fa_t& a, const float    b ) { return a * Vec3fa_t(b); }
  __forceinline const Vec3fa_t operator *( const float    a, const Vec3fa_t& b ) { return Vec3fa_t(a) * b; }
  __forceinline const Vec3fa_t operator /( const Vec3fa_t& a, const Vec3fa_t& b ) { return _mm512_div_ps(a.m512,b.m512); }
  __forceinline const Vec3fa_t operator /( const Vec3fa_t& a, const float    b ) { return _mm512_div_ps(a.m512,_mm512_set1_ps(b)); }
  __forceinline const Vec3fa_t operator /( const    float a, const Vec3fa_t& b ) { return _mm512_div_ps(_mm512_set1_ps(a),b.m512); }

  __forceinline const Vec3fa_t min( const Vec3fa_t& a, const Vec3fa_t& b ) { return _mm512_gmin_ps(a.m512,b.m512); }
  __forceinline const Vec3fa_t max( const Vec3fa_t& a, const Vec3fa_t& b ) { return _mm512_gmax_ps(a.m512,b.m512); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Ternary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline Vec3fa_t madd  ( const Vec3fa_t& a, const Vec3fa_t& b, const Vec3fa_t& c) { return _mm512_fmadd_ps(a,b,c); }
  __forceinline Vec3fa_t msub  ( const Vec3fa_t& a, const Vec3fa_t& b, const Vec3fa_t& c) { return _mm512_fmsub_ps(a,b,c); }
  __forceinline Vec3fa_t nmadd ( const Vec3fa_t& a, const Vec3fa_t& b, const Vec3fa_t& c) { return _mm512_fnmadd_ps(a,b,c); }
  __forceinline Vec3fa_t nmsub ( const Vec3fa_t& a, const Vec3fa_t& b, const Vec3fa_t& c) { return _mm512_fnmsub_ps(a,b,c); }


  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline Vec3fa& operator +=( Vec3fa& a, const Vec3fa_t& b ) { return a = a + b; }
  __forceinline Vec3fa& operator -=( Vec3fa& a, const Vec3fa_t& b ) { return a = a - b; }
  __forceinline Vec3fa& operator *=( Vec3fa& a, const float    b ) { return a = a * b; }
  __forceinline Vec3fa& operator /=( Vec3fa& a, const float    b ) { return a = a / b; }

  __forceinline Vec3fa_t& operator +=( Vec3fa_t& a, const Vec3fa_t& b ) { return a = a + b; }
  __forceinline Vec3fa_t& operator *=( Vec3fa_t& a, const float    b  ) { return a = a * b; }
  __forceinline Vec3fa_t& operator /=( Vec3fa_t& a, const float    b  ) { return a = a / b; }

  ////////////////////////////////////////////////////////////////////////////////
  /// Reductions
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline float reduce_add(const Vec3fa_t& a) { 
    __m512 b = _mm512_mask_blend_ps(0x8,a.m512,_mm512_set1_ps(0.0f));
    __m512 c = _mm512_add_ps(b, _mm512_swizzle_ps(b, _MM_SWIZ_REG_BADC));
    __m512 d = _mm512_add_ps(c, _mm512_swizzle_ps(c, _MM_SWIZ_REG_CDAB));
    return _mm512_cvtss_f32(d);
  } 

  __forceinline float reduce_mul(const Vec3fa_t& a) { 
    __m512 b = _mm512_mask_blend_ps(0x8,a.m512,_mm512_set1_ps(1.0f));
    __m512 c = _mm512_mul_ps(b, _mm512_swizzle_ps(b, _MM_SWIZ_REG_BADC));
    __m512 d = _mm512_mul_ps(c, _mm512_swizzle_ps(c, _MM_SWIZ_REG_CDAB));
    return _mm512_cvtss_f32(d);
  } 

  __forceinline float reduce_min(const Vec3fa_t& a) { 
    __m512 b = _mm512_mask_blend_ps(0x8,a.m512,_mm512_set1_ps(pos_inf));
    __m512 c = _mm512_gmin_ps(b, _mm512_swizzle_ps(b, _MM_SWIZ_REG_BADC));
    __m512 d = _mm512_gmin_ps(c, _mm512_swizzle_ps(c, _MM_SWIZ_REG_CDAB));
    return _mm512_cvtss_f32(d);
  } 

  __forceinline float reduce_max(const Vec3fa_t& a) { 
    __m512 b = _mm512_mask_blend_ps(0x8,a.m512,_mm512_set1_ps(neg_inf));
    __m512 c = _mm512_gmax_ps(b, _mm512_swizzle_ps(b, _MM_SWIZ_REG_BADC));
    __m512 d = _mm512_gmax_ps(c, _mm512_swizzle_ps(c, _MM_SWIZ_REG_CDAB));
    return _mm512_cvtss_f32(d);
  } 

  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline bool operator ==( const Vec3fa_t& a, const Vec3fa_t& b ) { return (_mm512_cmp_ps_mask(a,b,_MM_CMPINT_EQ) & 7) == 7; }
  __forceinline bool operator !=( const Vec3fa_t& a, const Vec3fa_t& b ) { return (_mm512_cmp_ps_mask(a,b,_MM_CMPINT_NE) & 7) != 0; }

  __forceinline Vec3ba eq_mask( const Vec3fa_t& a, const Vec3fa_t& b ) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_EQ); }
  __forceinline Vec3ba neq_mask(const Vec3fa_t& a, const Vec3fa_t& b ) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_NE); }
  __forceinline Vec3ba lt_mask( const Vec3fa_t& a, const Vec3fa_t& b ) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_LT); }
  __forceinline Vec3ba le_mask( const Vec3fa_t& a, const Vec3fa_t& b ) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_LE); }
  __forceinline Vec3ba gt_mask( const Vec3fa_t& a, const Vec3fa_t& b ) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_GT); }
  __forceinline Vec3ba ge_mask( const Vec3fa_t& a, const Vec3fa_t& b ) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_GE); }

 ////////////////////////////////////////////////////////////////////////////////
  /// Euclidian Space Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline float dot ( const Vec3fa_t& a, const Vec3fa_t& b ) {
    return reduce_add(a*b);
  }

  __forceinline Vec3fa_t cross ( const Vec3fa_t& a, const Vec3fa_t& b ) 
  {
    __m512 c = _mm512_mul_ps(b,_mm512_swizzle_ps(a,_MM_SWIZ_REG_DACB));
    c = _mm512_fmsub_ps(_mm512_swizzle_ps(b,_MM_SWIZ_REG_DACB),a,c);
    c = _mm512_swizzle_ps(c,_MM_SWIZ_REG_DACB);
    return c;
  }

  __forceinline float   length   ( const Vec3fa_t& a )                   { return sqrt(dot(a,a)); }
  __forceinline Vec3fa_t normalize( const Vec3fa_t& a )                   { return a*rsqrt(dot(a,a)); }
  __forceinline float   distance ( const Vec3fa_t& a, const Vec3fa_t& b ) { return length(a-b); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Select
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const Vec3fa_t select( bool s, const Vec3fa_t& t, const Vec3fa_t& f ) {
    return _mm512_mask_blend_ps(s ? 0xF : 0x0, f, t);
  }

  __forceinline const Vec3fa_t select( const Vec3ba& s, const Vec3fa_t& t, const Vec3fa_t& f ) {
    return _mm512_mask_blend_ps(s, f, t);
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Output Operators
  ////////////////////////////////////////////////////////////////////////////////

  inline std::ostream& operator<<(std::ostream& cout, const Vec3fa& a) {
    return cout << "(" << a.x << ", " << a.y << ", " << a.z << ", " << a.w << ")";
  }

  inline std::ostream& operator<<(std::ostream& cout, const Vec3fa_t& a) {
    const float * const f = (float*)&a.m512;
    return cout << "(" << f[0] << ", " << f[1] << ", " << f[2] << ", " << f[3] << ")";

  }

}

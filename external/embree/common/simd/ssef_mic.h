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

#ifndef __EMBREE_SSEF_MIC_H__
#define __EMBREE_SSEF_MIC_H__

namespace embree
{
  /* memory representation as 4 aligned floats */
  struct __aligned(16) ssef_m 
  {
    typedef sseb_m Mask;
    typedef ssei_m Int;
    typedef ssef_m Float;
    enum { size = 4 };  // number of SIMD elements
    union {
      float f[4]; 
      int i[4];
    };
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline ssef_m ( ) {}
    __forceinline ssef_m ( float a                           ) { f[0] = a;  f[1] = a; f[2] = a; f[3] = a; }
    __forceinline ssef_m ( float a, float b, float c, float d) { f[0] = a;  f[1] = b; f[2] = c; f[3] = d; }

    __forceinline ssef_m           ( const ssef_t& other ); 
    __forceinline ssef_m& operator=( const ssef_t& other );

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline ssef_m( ZeroTy   ) { f[0] = f[1] = f[2] = f[3] = zero; }
    __forceinline ssef_m( OneTy    ) { f[0] = f[1] = f[2] = f[3] = one; }
    __forceinline ssef_m( NegInfTy ) { f[0] = f[1] = f[2] = f[3] = neg_inf; }
    __forceinline ssef_m( PosInfTy ) { f[0] = f[1] = f[2] = f[3] = pos_inf; }
    __forceinline ssef_m( StepTy   ) { f[0] = 0.0f; f[1] = 1.0f; f[2] = 2.0f; f[3] = 3.0f; }
    __forceinline ssef_m( NaNTy    ) { f[0] = f[1] = f[2] = f[3] = nan; }

    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline const float& operator []( const size_t index ) const { assert(index < 4); return f[index]; }
    __forceinline       float& operator []( const size_t index )       { assert(index < 4); return f[index]; }
  };
  
  /*! 4-wide SSE float type emulated with 16-wide vectors. */
  struct __aligned(64) ssef_t 
  {
    union {
      __m512 m512; 
      float  f[4];
      int    i[4];
    };
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline ssef_t( const __m512 a ) : m512(a) {}
    __forceinline ssef_t           ( const ssef_t& other ) { m512 = other.m512; }
    __forceinline ssef_t& operator=( const ssef_t& other ) { m512 = other.m512; return *this; }

    __forceinline operator const __m512&( void ) const { return m512; }
    __forceinline operator       __m512&( void )       { return m512; }

    __forceinline ssef_t( const ssei_t& v );

    
  public:  
    __forceinline explicit ssef_t ( float a ) 
      : m512(_mm512_set1_ps(a)) {}

    __forceinline ssef_t ( const ssef_m& other ) { 
      m512 = _mm512_extload_ps(&other,_MM_UPCONV_PS_NONE,_MM_BROADCAST_4X16,_MM_HINT_NONE); 
    }
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline const float& operator []( const size_t index ) const { assert(index < 4); return f[index]; }
    __forceinline       float& operator []( const size_t index )       { assert(index < 4); return f[index]; }
  };

  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////
  
  /* __forceinline ssei_m::ssei_m ( const ssef_t& other ) { */
  /*   __m512i i = _mm512_cvtfxpnt_round_adjustps_epi32(other,_MM_ROUND_MODE_NEAREST,_MM_EXPADJ_NONE); */
  /*   _mm512_mask_extpackstorelo_epi32(this,0xf,i,_MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);  */
  /* } */

  __forceinline ssef_t::ssef_t( const ssei_t& v ) {
    m512 = _mm512_cvtfxpnt_round_adjustepi32_ps(v,_MM_ROUND_MODE_NEAREST,_MM_EXPADJ_NONE);
  }

  __forceinline ssei_t::ssei_t ( const ssef_t& other )
    {
      m512 = _mm512_cvtfxpnt_round_adjustps_epi32(other.m512,_MM_ROUND_NEAREST,_MM_EXPADJ_NONE);      
    }

  __forceinline ssef_m::ssef_m ( const ssef_t& other ) { 
    _mm512_mask_extpackstorelo_ps(this,0xf,other.m512,_MM_DOWNCONV_PS_NONE, _MM_HINT_NONE); 
  }
  
  __forceinline ssef_m& ssef_m::operator=( const ssef_t& other ) { 
    _mm512_mask_extpackstorelo_ps(this,0xf,other.m512,_MM_DOWNCONV_PS_NONE, _MM_HINT_NONE); 
    return *this; 
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Unary Operators
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline const ssef_t operator +( const ssef_t& a ) { return a; }
  __forceinline const ssef_t operator -( const ssef_t& a ) { return _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(a.m512), _mm512_set1_epi32(0x80000000))); }
  __forceinline const ssef_t abs       ( const ssef_t& a ) { return _mm512_gmaxabs_ps(a.m512,a.m512); } 
  __forceinline const ssef_t sign      ( const ssef_t& a ) { return _mm512_mask_blend_ps(_mm512_cmp_ps_mask(a,ssef_t(ssef_m(zero)),_MM_CMPINT_LT), ssef_t(ssef_m(one)), -ssef_t(ssef_m(one))); }
  __forceinline const ssef_t signmsk   ( const ssef_t& a ) { return _mm512_castsi512_ps(_mm512_and_epi32(_mm512_castps_si512(a.m512),_mm512_set1_epi32(0x80000000))); }
  
  __forceinline const ssef_t rcp  ( const ssef_t& a ) { return _mm512_rcp23_ps(a.m512); }
  __forceinline const ssef_t sqr  ( const ssef_t& a ) { return _mm512_mul_ps(a.m512,a.m512); }
  __forceinline const ssef_t sqrt ( const ssef_t& a ) { return _mm512_sqrt_ps(a.m512); }
  __forceinline const ssef_t rsqrt( const ssef_t& a ) { return _mm512_invsqrt_ps(a.m512); }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Binary Operators
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline const ssef_t operator +( const ssef_t& a, const ssef_t& b ) { return _mm512_add_ps(a.m512, b.m512); }
  __forceinline const ssef_t operator +( const ssef_t& a, const float & b ) { return a + ssef_t(b); }
  __forceinline const ssef_t operator +( const float & a, const ssef_t& b ) { return ssef_t(a) + b; }
  
  __forceinline const ssef_t operator -( const ssef_t& a, const ssef_t& b ) { return _mm512_sub_ps(a.m512, b.m512); }
  __forceinline const ssef_t operator -( const ssef_t& a, const float & b ) { return a - ssef_t(b); }
  __forceinline const ssef_t operator -( const float & a, const ssef_t& b ) { return ssef_t(a) - b; }
  
  __forceinline const ssef_t operator *( const ssef_t& a, const ssef_t& b ) { return _mm512_mul_ps(a.m512, b.m512); }
  __forceinline const ssef_t operator *( const ssef_t& a, const float & b ) { return a * ssef_t(b); }
  __forceinline const ssef_t operator *( const float & a, const ssef_t& b ) { return ssef_t(a) * b; }
  
  __forceinline const ssef_t operator /( const ssef_t& a, const ssef_t& b ) { return _mm512_div_ps(a.m512,b.m512); }
  __forceinline const ssef_t operator /( const ssef_t& a, const float & b ) { return a/ssef_t(b); }
  __forceinline const ssef_t operator /( const float & a, const ssef_t& b ) { return ssef_t(a)/b; }
  
  __forceinline const ssef_t operator^( const ssef_t& a, const ssef_t& b ) { return _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(a.m512),_mm512_castps_si512(b.m512))); }
  __forceinline const ssef_t operator^( const ssef_t& a, const ssei_t& b ) { return _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(a.m512),b.m512)); }
  
  __forceinline const ssef_t min( const ssef_t& a, const ssef_t& b ) { return _mm512_gmin_ps(a.m512,b.m512); }
  __forceinline const ssef_t min( const ssef_t& a, const float & b ) { return _mm512_gmin_ps(a.m512,ssef_t(b)); }
  __forceinline const ssef_t min( const float & a, const ssef_t& b ) { return _mm512_gmin_ps(ssef_t(a),b.m512); }
  
  __forceinline const ssef_t max( const ssef_t& a, const ssef_t& b ) { return _mm512_gmax_ps(a.m512,b.m512); }
  __forceinline const ssef_t max( const ssef_t& a, const float & b ) { return _mm512_gmax_ps(a.m512,ssef_t(b)); }
  __forceinline const ssef_t max( const float & a, const ssef_t& b ) { return _mm512_gmax_ps(ssef_t(a),b.m512); }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Ternary Operators
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline const ssef_t madd  ( const ssef_t& a, const ssef_t& b, const ssef_t& c) { return _mm512_fmadd_ps(a,b,c); }
  __forceinline const ssef_t msub  ( const ssef_t& a, const ssef_t& b, const ssef_t& c) { return _mm512_fmsub_ps(a,b,c); }
  __forceinline const ssef_t nmadd ( const ssef_t& a, const ssef_t& b, const ssef_t& c) { return _mm512_fnmadd_ps(a,b,c); }
  __forceinline const ssef_t nmsub ( const ssef_t& a, const ssef_t& b, const ssef_t& c) { return _mm512_fnmsub_ps(a,b,c); }

  ////////////////////////////////////////////////////////////////////////////////
  /// Assignment Operators
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline ssef_t& operator +=( ssef_t& a, const ssef_t& b ) { return a = a + b; }
  __forceinline ssef_t& operator +=( ssef_t& a, const float & b ) { return a = a + b; }
  
  __forceinline ssef_t& operator -=( ssef_t& a, const ssef_t& b ) { return a = a - b; }
  __forceinline ssef_t& operator -=( ssef_t& a, const float & b ) { return a = a - b; }
  
  __forceinline ssef_t& operator *=( ssef_t& a, const ssef_t& b ) { return a = a * b; }
  __forceinline ssef_t& operator *=( ssef_t& a, const float & b ) { return a = a * b; }
  
  __forceinline ssef_t& operator /=( ssef_t& a, const ssef_t& b ) { return a = a / b; }
  __forceinline ssef_t& operator /=( ssef_t& a, const float & b ) { return a = a / b; }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators + Select
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline const sseb_t operator ==( const ssef_t& a, const ssef_t& b ) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_EQ); }
  __forceinline const sseb_t operator ==( const ssef_t& a, const float & b ) { return a == ssef_t(b); }
  __forceinline const sseb_t operator ==( const float & a, const ssef_t& b ) { return ssef_t(a) == b; }
  
  __forceinline const sseb_t operator !=( const ssef_t& a, const ssef_t& b ) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_NE); }
  __forceinline const sseb_t operator !=( const ssef_t& a, const float & b ) { return a != ssef_t(b); }
  __forceinline const sseb_t operator !=( const float & a, const ssef_t& b ) { return ssef_t(a) != b; }
  
  __forceinline const sseb_t operator < ( const ssef_t& a, const ssef_t& b ) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_LT); }
  __forceinline const sseb_t operator < ( const ssef_t& a, const float & b ) { return a <  ssef_t(b); }
  __forceinline const sseb_t operator < ( const float & a, const ssef_t& b ) { return ssef_t(a) <  b; }
  
  __forceinline const sseb_t operator >=( const ssef_t& a, const ssef_t& b ) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_GE); }
  __forceinline const sseb_t operator >=( const ssef_t& a, const float & b ) { return a >= ssef_t(b); }
  __forceinline const sseb_t operator >=( const float & a, const ssef_t& b ) { return ssef_t(a) >= b; }
  
  __forceinline const sseb_t operator > ( const ssef_t& a, const ssef_t& b ) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_GT); }
  __forceinline const sseb_t operator > ( const ssef_t& a, const float & b ) { return a >  ssef_t(b); }
  __forceinline const sseb_t operator > ( const float & a, const ssef_t& b ) { return ssef_t(a) >  b; }
  
  __forceinline const sseb_t operator <=( const ssef_t& a, const ssef_t& b ) { return _mm512_cmp_ps_mask(a,b,_MM_CMPINT_LE); }
  __forceinline const sseb_t operator <=( const ssef_t& a, const float & b ) { return a <= ssef_t(b); }
  __forceinline const sseb_t operator <=( const float & a, const ssef_t& b ) { return ssef_t(a) <= b; }
  
  __forceinline const ssef_t select( const sseb_t& mask, const ssef_t& t, const ssef_t& f ) { return _mm512_mask_blend_ps(mask.v, f, t); }


  
  ////////////////////////////////////////////////////////////////////////////////
  /// Rounding Functions
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline const ssef_t round_even( const ssef_t& a ) { return _mm512_round_ps(a, _MM_ROUND_MODE_NEAREST    , _MM_EXPADJ_NONE); }
  __forceinline const ssef_t round_down( const ssef_t& a ) { return _mm512_round_ps(a, _MM_ROUND_MODE_DOWN       , _MM_EXPADJ_NONE); }
  __forceinline const ssef_t round_up  ( const ssef_t& a ) { return _mm512_round_ps(a, _MM_ROUND_MODE_UP         , _MM_EXPADJ_NONE); }
  __forceinline const ssef_t round_zero( const ssef_t& a ) { return _mm512_round_ps(a, _MM_ROUND_MODE_TOWARD_ZERO, _MM_EXPADJ_NONE); }
  __forceinline const ssef_t floor     ( const ssef_t& a ) { return _mm512_round_ps(a, _MM_ROUND_MODE_DOWN       , _MM_EXPADJ_NONE); }
  __forceinline const ssef_t ceil      ( const ssef_t& a ) { return _mm512_round_ps(a, _MM_ROUND_MODE_UP         , _MM_EXPADJ_NONE); }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Movement/Shifting/Shuffling Functions
  ////////////////////////////////////////////////////////////////////////////////
  
  //__forceinline ssef_t unpacklo( const ssef_t& a, const ssef_t& b ) { return _mm512_unpacklo_ps(a.m512, b.m512); }
  //__forceinline ssef_t unpackhi( const ssef_t& a, const ssef_t& b ) { return _mm512_unpackhi_ps(a.m512, b.m512); }
  
  template<size_t i0, size_t i1, size_t i2, size_t i3> __forceinline const ssef_t shuffle( const ssef_t& b ) {
    return _mm512_castsi512_ps(_mm512_shuffle_epi32(_mm512_castps_si512(b),(_MM_PERM_ENUM)_MM_SHUFFLE(i3, i2, i1, i0)));
  }
  
  //template<size_t i0, size_t i1, size_t i2, size_t i3> __forceinline const ssef_t shuffle( const ssef_t& a, const ssef_t& b ) {
  //  return _mm512_shuffle_ps(a, b, _MM512_SHUFFLE(i3, i2, i1, i0));
  //}
  
  //__forceinline const ssef_t shuffle8(const ssef_t& a, const ssei& shuf) { 
  //  return _mm512_castsi512_ps(_mm512_shuffle_epi8(_mm512_castps_si512(a), shuf)); 
  //}

  template<size_t i0> __forceinline const ssef_t shuffle( const ssef_t& b ) {
    return shuffle<i0,i0,i0,i0>(b);
  }

  template<size_t i> __forceinline float extract( const ssef_t& a ) { return a[i]; }
  //template<size_t dst, size_t src, size_t clr> __forceinline const ssef_t insert( const ssef_t& a, const ssef_t& b ) { return _mm512_insert_ps(a, b, (dst << 4) | (src << 6) | clr); }
  template<size_t dst, size_t src> __forceinline const ssef_t insert( const ssef_t& a, const ssef_t& b ) { ssef_t c = a; c[dst] = b[src]; return c; }
  template<size_t dst>             __forceinline const ssef_t insert( const ssef_t& a, const float b ) { ssef_t c = a; c[dst] = b; return c; }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Transpose
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline void transpose(const ssef_t& r0, const ssef_t& r1, const ssef_t& r2, const ssef_t& r3, ssef_m& c0, ssef_m& c1, ssef_m& c2, ssef_m& c3)
  {
    c0[0] = r0[0]; c0[1] = r1[0]; c0[2] = r2[0]; c0[3] = r3[0]; 
    c1[0] = r0[1]; c1[1] = r1[1]; c1[2] = r2[1]; c1[3] = r3[1]; 
    c2[0] = r0[2]; c2[1] = r1[2]; c2[2] = r2[2]; c2[3] = r3[2]; 
    c3[0] = r0[3]; c3[1] = r1[3]; c3[2] = r2[3]; c3[3] = r3[3]; 
  }
  
  __forceinline void transpose(const ssef_t& r0, const ssef_t& r1, const ssef_t& r2, const ssef_t& r3, ssef_m& c0, ssef_m& c1, ssef_m& c2) 
  {
    c0[0] = r0[0]; c0[1] = r1[0]; c0[2] = r2[0]; c0[3] = r3[0]; 
    c1[0] = r0[1]; c1[1] = r1[1]; c1[2] = r2[1]; c1[3] = r3[1]; 
    c2[0] = r0[2]; c2[1] = r1[2]; c2[2] = r2[2]; c2[3] = r3[2]; 
  }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Reductions
  ////////////////////////////////////////////////////////////////////////////////
  
  __forceinline const ssef_t vreduce_min(const ssef_t& v) { ssef_t h = min(shuffle<1,0,3,2>(v),v); return min(shuffle<2,3,0,1>(h),h); }
  __forceinline const ssef_t vreduce_max(const ssef_t& v) { ssef_t h = max(shuffle<1,0,3,2>(v),v); return max(shuffle<2,3,0,1>(h),h); }
  __forceinline const ssef_t vreduce_add(const ssef_t& v) { ssef_t h = shuffle<1,0,3,2>(v)   + v ; return shuffle<2,3,0,1>(h)   + h ; }
  
  __forceinline float reduce_min(const ssef_t& v) { return _mm512_cvtss_f32(vreduce_min(v)); }
  __forceinline float reduce_max(const ssef_t& v) { return _mm512_cvtss_f32(vreduce_max(v)); }
  __forceinline float reduce_add(const ssef_t& v) { return _mm512_cvtss_f32(vreduce_add(v)); }
  
  __forceinline size_t select_min(const ssef_t& v) { return bitscan(movemask(v == vreduce_min(v))); }
  __forceinline size_t select_max(const ssef_t& v) { return bitscan(movemask(v == vreduce_max(v))); }
  
  __forceinline size_t select_min(const sseb_t& valid, const ssef_t& v) { 
    const ssef_t a = select(valid,v,ssef_m(pos_inf)); 
    return bitscan(movemask(valid & (a == vreduce_min(a)))); 
  }
  __forceinline size_t select_max(const sseb_t& valid, const ssef_t& v) { 
    const ssef_t a = select(valid,v,ssef_m(neg_inf));
    return bitscan(movemask(valid & (a == vreduce_max(a)))); 
  }

 ////////////////////////////////////////////////////////////////////////////////
  /// Memory load and store operations
  ////////////////////////////////////////////////////////////////////////////////

  __forceinline ssef_t load4f( const void* const a ) {
    return _mm512_extload_ps(a,_MM_UPCONV_PS_NONE,_MM_BROADCAST_4X16,_MM_HINT_NONE);
  }

  __forceinline void store4f( void* ptr, const ssef_t& v ) {
    assert((size_t)ptr % 16 == 0); 
    _mm512_mask_extpackstorelo_ps(ptr,0xf,v.m512,_MM_DOWNCONV_PS_NONE, _MM_HINT_NONE); 
  }

  __forceinline void storeu4f( void* ptr, const ssef_t& v ) {
    _mm512_mask_extpackstorelo_ps(ptr,0xf,v.m512,_MM_DOWNCONV_PS_NONE, _MM_HINT_NONE);
    _mm512_mask_extpackstorehi_ps(ptr,0xf,v.m512,_MM_DOWNCONV_PS_NONE, _MM_HINT_NONE);
  }

  __forceinline void store4f( const sseb_t& mask, void* ptr, const ssef_t& f ) { 
    store4f(ptr,select(mask,f,*(ssef_m*)ptr));
  }

  __forceinline ssef_t load4f_nt(void* ptr) {
    return _mm512_extload_ps(ptr,_MM_UPCONV_PS_NONE,_MM_BROADCAST_4X16,_MM_HINT_NT);
  }

  __forceinline void store4f_nt(void* ptr, const ssef_t& v) {
    assert((size_t)ptr % 16 == 0); 
    _mm512_mask_extpackstorelo_ps(ptr,0xf,v.m512,_MM_DOWNCONV_PS_NONE, _MM_HINT_NT); 
  }
  
  ////////////////////////////////////////////////////////////////////////////////
  /// Output Operators
  ////////////////////////////////////////////////////////////////////////////////
  
  inline std::ostream& operator<<(std::ostream& cout, const ssef_m& a) {
    return cout << "<" << a[0] << ", " << a[1] << ", " << a[2] << ", " << a[3] << ">";
  }  
}

#endif

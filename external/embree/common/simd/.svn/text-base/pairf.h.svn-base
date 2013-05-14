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

#ifndef __EMBREE_PAIRF_H__
#define __EMBREE_PAIRF_H__

namespace embree
{
  /*! Doubling SIMD float types. */
  template<typename T> struct pairf
  {
    typedef pairb<typename T::Mask> Mask;  // mask type for us
    typedef pairi<typename T::Int > Int ;  // int type for us
    enum   { size = 2*T::size };           // number of SIMD elements
    T l,h;                                 // low and high data element

    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline pairf           ( ) {}
    __forceinline pairf           ( const pairf& a ) { l = a.l; h = a.h; }
    __forceinline pairf& operator=( const pairf& a ) { l = a.l; h = a.h; return *this; }

    __forceinline pairf           ( const T& a             ) : l(a), h(a) {}
    __forceinline pairf           ( const T& a, const T& b ) : l(a), h(b) {}

    __forceinline explicit pairf  ( const float* const a ) : l(a), h(a+T::size) {}
    __forceinline pairf           ( float  a ) : l(a), h(a) {}
    __forceinline pairf           ( float  a, float  b) : l(a), h(b) {}
    __forceinline pairf           ( float  a, float  b, float  c, float  d) : l(a,b), h(c,d) {}
    __forceinline pairf           ( float  a, float  b, float  c, float  d, float  e, float  f, float  g, float  h) : l(a,b,c,d), h(e,f,g,h) {}

    __forceinline explicit pairf  ( const Int& a ) : l(a.l), h(a.h) {}

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline pairf( ZeroTy   ) : l(zero   ), h(zero   ) {}
    __forceinline pairf( OneTy    ) : l(one    ), h(one    ) {}
    __forceinline pairf( PosInfTy ) : l(pos_inf), h(pos_inf) {}
    __forceinline pairf( NegInfTy ) : l(neg_inf), h(neg_inf) {}
    __forceinline pairf( StepTy )   : l(step   ), h(step   ) { h += T::size; }
    __forceinline pairf( NaNTy    ) : l(nan    ), h(nan    ) {}

    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline const float& operator []( const size_t i ) const { assert(i < size); return ((float*)this)[i]; }
    __forceinline       float& operator []( const size_t i )       { assert(i < size); return ((float*)this)[i]; }

    ////////////////////////////////////////////////////////////////////////////////
    /// Unary Operators
    ////////////////////////////////////////////////////////////////////////////////
    
    friend __forceinline const pairf operator +( const pairf& a ) { return a; }
    friend __forceinline const pairf operator -( const pairf& a ) { return pairf(       -a.l ,        -a.h ); }
    friend __forceinline const pairf abs       ( const pairf& a ) { return pairf(    abs(a.l),     abs(a.h)); }
    friend __forceinline const pairf sign      ( const pairf& a ) { return pairf(   sign(a.l)   , sign(a.h)); }
    friend __forceinline const pairf signmsk   ( const pairf& a ) { return pairf(signmsk(a.l), signmsk(a.h)); }
    
    friend __forceinline const pairf rcp       ( const pairf& a ) { return pairf(  rcp(a.l),  rcp(a.h)); }
    friend __forceinline const pairf sqr       ( const pairf& a ) { return pairf(  sqr(a.l),  sqr(a.h)); }
    friend __forceinline const pairf sqrt      ( const pairf& a ) { return pairf( sqrt(a.l), sqrt(a.h)); }
    friend __forceinline const pairf rsqrt     ( const pairf& a ) { return pairf(rsqrt(a.l),rsqrt(a.h)); }
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Binary Operators
    ////////////////////////////////////////////////////////////////////////////////
    
    friend __forceinline const pairf operator +( const pairf& a, const pairf& b ) { return pairf(a.l + b.l, a.h + b.h); }
    friend __forceinline const pairf operator +( const pairf& a, const float  b ) { return a + pairf(b); }
    friend __forceinline const pairf operator +( const float  a, const pairf& b ) { return pairf(a) + b; }
    
    friend __forceinline const pairf operator -( const pairf& a, const pairf& b ) { return pairf(a.l - b.l, a.h - b.h); }
    friend __forceinline const pairf operator -( const pairf& a, const float  b ) { return a - pairf(b); }
    friend __forceinline const pairf operator -( const float  a, const pairf& b ) { return pairf(a) - b; }
    
    friend __forceinline const pairf operator *( const pairf& a, const pairf& b ) { return pairf(a.l * b.l, a.h * b.h); }
    friend __forceinline const pairf operator *( const pairf& a, const float  b ) { return a * pairf(b); }
    friend __forceinline const pairf operator *( const float  a, const pairf& b ) { return pairf(a) * b; }
    
    friend __forceinline const pairf operator /( const pairf& a, const pairf& b ) { return a * rcp(b); }
    friend __forceinline const pairf operator /( const pairf& a, const float  b ) { return a * rcp(b); }
    friend __forceinline const pairf operator /( const float  a, const pairf& b ) { return a * rcp(b); }
    
    friend __forceinline const pairf operator^( const pairf& a, const pairf& b ) { return pairf(a.l ^ b.l, a.h ^ b.h); }
    friend __forceinline const pairf operator^( const pairf& a, const Int  & b ) { return pairf(a.l ^ b.l, a.h ^ b.h); }
    
    friend __forceinline const pairf min( const pairf& a, const pairf& b ) { return pairf(min(a.l,b.l), min(a.h,b.h)); }
    friend __forceinline const pairf min( const pairf& a, const float  b ) { return min(a,pairf(b)); }
    friend __forceinline const pairf min( const float  a, const pairf& b ) { return min(pairf(a),b); }
    
    friend __forceinline const pairf max( const pairf& a, const pairf& b ) { return pairf(max(a.l,b.l), max(a.h,b.h)); }
    friend __forceinline const pairf max( const pairf& a, const float  b ) { return max(a,pairf(b)); }
    friend __forceinline const pairf max( const float  a, const pairf& b ) { return max(pairf(a),b); }
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Assignment Operators
    ////////////////////////////////////////////////////////////////////////////////
    
    friend __forceinline pairf& operator +=( pairf& a, const pairf& b ) { return a = a + b; }
    friend __forceinline pairf& operator +=( pairf& a, const float  b ) { return a = a + b; }
    
    friend __forceinline pairf& operator -=( pairf& a, const pairf& b ) { return a = a - b; }
    friend __forceinline pairf& operator -=( pairf& a, const float  b ) { return a = a - b; }
    
    friend __forceinline pairf& operator *=( pairf& a, const pairf& b ) { return a = a * b; }
    friend __forceinline pairf& operator *=( pairf& a, const float  b ) { return a = a * b; }
    
    friend __forceinline pairf& operator /=( pairf& a, const pairf& b ) { return a = a / b; }
    friend __forceinline pairf& operator /=( pairf& a, const float  b ) { return a = a / b; }
    
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Comparison Operators + Select
    ////////////////////////////////////////////////////////////////////////////////
    
    friend __forceinline const Mask operator ==( const pairf& a, const pairf& b ) { return Mask(a.l == b.l, a.h == b.h); }
    friend __forceinline const Mask operator ==( const pairf& a, const float  b ) { return a == pairf(b); }
    friend __forceinline const Mask operator ==( const float  a, const pairf& b ) { return pairf(a) == b; }
    
    friend __forceinline const Mask operator !=( const pairf& a, const pairf& b ) { return Mask(a.l != b.l, a.h != b.h); }
    friend __forceinline const Mask operator !=( const pairf& a, const float  b ) { return a != pairf(b); }
    friend __forceinline const Mask operator !=( const float  a, const pairf& b ) { return pairf(a) != b; }
    
    friend __forceinline const Mask operator < ( const pairf& a, const pairf& b ) { return Mask(a.l <  b.l, a.h <  b.h); }
    friend __forceinline const Mask operator < ( const pairf& a, const float  b ) { return a <  pairf(b); }
    friend __forceinline const Mask operator < ( const float  a, const pairf& b ) { return pairf(a) <  b; }
    
    friend __forceinline const Mask operator >=( const pairf& a, const pairf& b ) { return Mask(a.l >= b.l, a.h >= b.h); }
    friend __forceinline const Mask operator >=( const pairf& a, const float  b ) { return a >= pairf(b); }
    friend __forceinline const Mask operator >=( const float  a, const pairf& b ) { return pairf(a) >= b; }
    
    friend __forceinline const Mask operator > ( const pairf& a, const pairf& b ) { return Mask(a.l >  b.l, a.h >  b.h); }
    friend __forceinline const Mask operator > ( const pairf& a, const float  b ) { return a >  pairf(b); }
    friend __forceinline const Mask operator > ( const float  a, const pairf& b ) { return pairf(a) >  b; }
    
    friend __forceinline const Mask operator <=( const pairf& a, const pairf& b ) { return Mask(a.l <= b.l, a.h <= b.h); }
    friend __forceinline const Mask operator <=( const pairf& a, const float  b ) { return a <= pairf(b); }
    friend __forceinline const Mask operator <=( const float  a, const pairf& b ) { return pairf(a) <= b; }
    
    friend __forceinline const pairf select ( const Mask& m, const pairf& a, const pairf& b ) { 
      return pairf(select(m.l,a.l,b.l),select(m.h,a.h,b.h)); 
    }
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Rounding Functions
    ////////////////////////////////////////////////////////////////////////////////
    
    friend __forceinline const pairf round_even( const pairf& a ) { return pairf(round_even(a.l), round_even(a.h)); }
    friend __forceinline const pairf round_down( const pairf& a ) { return pairf(round_down(a.l), round_down(a.h)); }
    friend __forceinline const pairf round_up  ( const pairf& a ) { return pairf(round_up  (a.l), round_up  (a.h)); }
    friend __forceinline const pairf round_zero( const pairf& a ) { return pairf(round_zero(a.l), round_zero(a.h)); }
    friend __forceinline const pairf floor     ( const pairf& a ) { return pairf(floor     (a.l), floor     (a.h)); }
    friend __forceinline const pairf ceil      ( const pairf& a ) { return pairf(ceil      (a.l), ceil      (a.h)); }
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Reductions
    ////////////////////////////////////////////////////////////////////////////////
    
    friend __forceinline const pairf vreduce_min(const pairf& v) { const T l = vreduce_min(v.l); const T h = vreduce_min(v.h); return pairf(min(l,h)); }
    friend __forceinline const pairf vreduce_max(const pairf& v) { const T l = vreduce_max(v.l); const T h = vreduce_max(v.h); return pairf(max(l,h)); }
    friend __forceinline const pairf vreduce_add(const pairf& v) { const T l = vreduce_add(v.l); const T h = vreduce_add(v.h); return pairf(l+h); }

    friend __forceinline float reduce_min(const pairf& v) { return min(reduce_min(v.l), reduce_min(v.h)); }
    friend __forceinline float reduce_max(const pairf& v) { return max(reduce_max(v.l), reduce_max(v.h)); }
    friend __forceinline float reduce_add(const pairf& v) { return     reduce_add(v.l)+ reduce_add(v.h) ; }

    friend __forceinline size_t select_min(const pairf& v) { return __bsf(movemask(v == vreduce_min(v))); }
    friend __forceinline size_t select_max(const pairf& v) { return __bsf(movemask(v == vreduce_max(v))); }

    friend __forceinline size_t select_min(const Mask& valid, const pairf& v) { 
      const pairf a = select(valid,v,pairf(pos_inf)); 
      return __bsf(movemask(valid & (a == vreduce_min(a)))); 
    }
    friend __forceinline size_t select_max(const Mask& valid, const pairf& v) { 
      const pairf a = select(valid,v,pairf(neg_inf)); 
      return __bsf(movemask(valid & (a == vreduce_max(a)))); 
    }

    friend __forceinline const T extract0(const pairf& v) { return v.l; }
    friend __forceinline const T extract1(const pairf& v) { return v.h; }
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Output Operators
    ////////////////////////////////////////////////////////////////////////////////
    
    friend inline std::ostream& operator<<(std::ostream& cout, const pairf& a) {
      return cout << "<" << a.l << ", " << a.h << ">";
    }
  };
}

#endif

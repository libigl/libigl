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

#ifndef __EMBREE_PAIRI_H__
#define __EMBREE_PAIRI_H__

namespace embree
{
  /*! Doubling SIMD int types. */
  template<typename T> struct pairi
  {
    typedef pairb<typename T::Mask> Mask;  // mask type for us
    enum   { size = 2*T::size };           // number of SIMD elements
    T l,h;                                 // low and high data element

    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline pairi           ( ) {}
    __forceinline pairi           ( const pairi& a ) { l = a.l; h = a.h; }
    __forceinline pairi& operator=( const pairi& a ) { l = a.l; h = a.h; return *this; }

    __forceinline pairi           ( const T& a             ) : l(a), h(a) {}
    __forceinline pairi           ( const T& a, const T& b ) : l(a), h(b) {}

    __forceinline explicit pairi  ( const int32* const a ) : l(a), h(a+T::size) {}
    __forceinline pairi           ( int32  a ) : l(a), h(a) {}
    __forceinline pairi           ( int32  a, int32  b) : l(a), h(b) {}
    __forceinline pairi           ( int32  a, int32  b, int32  c, int32  d) : l(a,b), h(c,d) {}
    __forceinline pairi           ( int32  a, int32  b, int32  c, int32  d, int32  e, int32  f, int32  g, int32  h) : l(a,b,c,d), h(e,f,g,h) {}
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////
    
    __forceinline pairi( ZeroTy   ) : l(zero   ), h(zero   ) {}
    __forceinline pairi( OneTy    ) : l(one    ), h(one    ) {}
    __forceinline pairi( PosInfTy ) : l(pos_inf), h(pos_inf) {}
    __forceinline pairi( NegInfTy ) : l(neg_inf), h(neg_inf) {}
    __forceinline pairi( StepTy )   : l(step   ), h(step   ) { h += T::size; }
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline const int32& operator []( const size_t i ) const { assert(i < size); return ((int32*)this)[i]; }
    __forceinline       int32& operator []( const size_t i )       { assert(i < size); return ((int32*)this)[i]; }

    ////////////////////////////////////////////////////////////////////////////////
    /// Unary Operators
    ////////////////////////////////////////////////////////////////////////////////
    
    friend __forceinline const pairi operator +( const pairi& a ) { return a; }
    friend __forceinline const pairi operator -( const pairi& a ) { return pairi(-a.l,-a.h); }
    friend __forceinline const pairi abs       ( const pairi& a ) { return pairi(abs(a.l),abs(a.h)); }
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Binary Operators
    ////////////////////////////////////////////////////////////////////////////////
    
    friend __forceinline const pairi operator +( const pairi& a, const pairi& b ) { return pairi(a.l + b.l, a.h + b.h); }
    friend __forceinline const pairi operator +( const pairi& a, const int32  b ) { return a + pairi(b); }
    friend __forceinline const pairi operator +( const int32  a, const pairi& b ) { return pairi(a) + b; }
    
    friend __forceinline const pairi operator -( const pairi& a, const pairi& b ) { return pairi(a.l - b.l, a.h - b.h); }
    friend __forceinline const pairi operator -( const pairi& a, const int32  b ) { return a - pairi(b); }
    friend __forceinline const pairi operator -( const int32  a, const pairi& b ) { return pairi(a) - b; }
    
    friend __forceinline const pairi operator *( const pairi& a, const pairi& b ) { return pairi(a.l * b.l, a.h * b.h); }
    friend __forceinline const pairi operator *( const pairi& a, const int32  b ) { return a * pairi(b); }
    friend __forceinline const pairi operator *( const int32  a, const pairi& b ) { return pairi(a) * b; }
    
    friend __forceinline const pairi operator &( const pairi& a, const pairi& b ) { return pairi(a.l & b.l, a.h & b.h); }
    friend __forceinline const pairi operator &( const pairi& a, const int32  b ) { return a & pairi(b); }
    friend __forceinline const pairi operator &( const int32  a, const pairi& b ) { return pairi(a) & b; }
    
    friend __forceinline const pairi operator |( const pairi& a, const pairi& b ) { return pairi(a.l | b.l, a.h | b.h); }
    friend __forceinline const pairi operator |( const pairi& a, const int32  b ) { return a | pairi(b); }
    friend __forceinline const pairi operator |( const int32  a, const pairi& b ) { return pairi(a) | b; }
    
    friend __forceinline const pairi operator ^( const pairi& a, const pairi& b ) { return pairi(a.l ^ b.l, a.h ^ b.h); }
    friend __forceinline const pairi operator ^( const pairi& a, const int32  b ) { return a ^ pairi(b); }
    friend __forceinline const pairi operator ^( const int32  a, const pairi& b ) { return pairi(a) ^ b; }
    
    friend __forceinline const pairi operator <<( const pairi& a, const int32  b ) { return pairi(a.l <<  b, a.h <<  b); }
    friend __forceinline const pairi operator >>( const pairi& a, const int32  b ) { return pairi(a.l >>  b, a.h >>  b); }
    
    friend __forceinline const pairi sra ( const pairi& a, const int32 b ) { return pairi(sra(a.l,b), sra(a.h,b)); }
    friend __forceinline const pairi srl ( const pairi& a, const int32 b ) { return pairi(srl(a.l,b), srl(a.h,b)); }
    
    friend __forceinline const pairi min( const pairi& a, const pairi& b ) { return pairi(min(a.l,b.l), min(a.h,b.h)); }
    friend __forceinline const pairi min( const pairi& a, const int32  b ) { return min(a,pairi(b)); }
    friend __forceinline const pairi min( const int32  a, const pairi& b ) { return min(pairi(a),b); }
    
    friend __forceinline const pairi max( const pairi& a, const pairi& b ) { return pairi(max(a.l,b.l), max(a.h,b.h)); }
    friend __forceinline const pairi max( const pairi& a, const int32  b ) { return max(a,pairi(b)); }
    friend __forceinline const pairi max( const int32  a, const pairi& b ) { return max(pairi(a),b); }
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Assignment Operators
    ////////////////////////////////////////////////////////////////////////////////
    
    friend __forceinline pairi& operator +=( pairi& a, const pairi& b ) { return a = a + b; }
    friend __forceinline pairi& operator +=( pairi& a, const int32  b ) { return a = a + b; }
    
    friend __forceinline pairi& operator -=( pairi& a, const pairi& b ) { return a = a - b; }
    friend __forceinline pairi& operator -=( pairi& a, const int32  b ) { return a = a - b; }
    
    friend __forceinline pairi& operator *=( pairi& a, const pairi& b ) { return a = a * b; }
    friend __forceinline pairi& operator *=( pairi& a, const int32  b ) { return a = a * b; }
    
    friend __forceinline pairi& operator &=( pairi& a, const pairi& b ) { return a = a & b; }
    friend __forceinline pairi& operator &=( pairi& a, const int32  b ) { return a = a & b; }
    
    friend __forceinline pairi& operator |=( pairi& a, const pairi& b ) { return a = a | b; }
    friend __forceinline pairi& operator |=( pairi& a, const int32  b ) { return a = a | b; }
    
    friend __forceinline pairi& operator <<=( pairi& a, const int32  b ) { return a = a << b; }
    friend __forceinline pairi& operator >>=( pairi& a, const int32  b ) { return a = a >> b; }
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Comparison Operators + Select
    ////////////////////////////////////////////////////////////////////////////////
    
    friend __forceinline const Mask operator ==( const pairi& a, const pairi& b ) { return Mask(a.l == b.l, a.h == b.h); }
    friend __forceinline const Mask operator ==( const pairi& a, const int32  b ) { return a == pairi(b); }
    friend __forceinline const Mask operator ==( const int32  a, const pairi& b ) { return pairi(a) == b; }
    
    friend __forceinline const Mask operator !=( const pairi& a, const pairi& b ) { return Mask(a.l != b.l, a.h != b.h); }
    friend __forceinline const Mask operator !=( const pairi& a, const int32  b ) { return a != pairi(b); }
    friend __forceinline const Mask operator !=( const int32  a, const pairi& b ) { return pairi(a) != b; }
    
    friend __forceinline const Mask operator < ( const pairi& a, const pairi& b ) { return Mask(a.l <  b.l, a.h <  b.h); }
    friend __forceinline const Mask operator < ( const pairi& a, const int32  b ) { return a <  pairi(b); }
    friend __forceinline const Mask operator < ( const int32  a, const pairi& b ) { return pairi(a) <  b; }
    
    friend __forceinline const Mask operator >=( const pairi& a, const pairi& b ) { return Mask(a.l >= b.l, a.h >= b.h); }
    friend __forceinline const Mask operator >=( const pairi& a, const int32  b ) { return a >= pairi(b); }
    friend __forceinline const Mask operator >=( const int32  a, const pairi& b ) { return pairi(a) >= b; }
    
    friend __forceinline const Mask operator > ( const pairi& a, const pairi& b ) { return Mask(a.l >  b.l, a.h >  b.h); }
    friend __forceinline const Mask operator > ( const pairi& a, const int32  b ) { return a >  pairi(b); }
    friend __forceinline const Mask operator > ( const int32  a, const pairi& b ) { return pairi(a) >  b; }
    
    friend __forceinline const Mask operator <=( const pairi& a, const pairi& b ) { return Mask(a.l <= b.l, a.h <= b.h); }
    friend __forceinline const Mask operator <=( const pairi& a, const int32  b ) { return a <= pairi(b); }
    friend __forceinline const Mask operator <=( const int32  a, const pairi& b ) { return pairi(a) <= b; }
    
    friend __forceinline const pairi select ( const Mask& m, const pairi& a, const pairi& b ) { 
      return pairi(select(m.l,a.l,b.l),select(m.h,a.h,b.h)); 
    }
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Reductions
    ////////////////////////////////////////////////////////////////////////////////
    
    friend __forceinline const pairi vreduce_min(const pairi& v) { const T l = vreduce_min(v.l); const T h = vreduce_min(v.h); return pairi(min(l,h)); }
    friend __forceinline const pairi vreduce_max(const pairi& v) { const T l = vreduce_max(v.l); const T h = vreduce_max(v.h); return pairi(max(l,h)); }
    friend __forceinline const pairi vreduce_add(const pairi& v) { const T l = vreduce_add(v.l); const T h = vreduce_add(v.h); return pairi(l+h); }

    friend __forceinline float reduce_min(const pairi& v) { return min(reduce_min(v.l), reduce_min(v.h)); }
    friend __forceinline float reduce_max(const pairi& v) { return max(reduce_max(v.l), reduce_max(v.h)); }
    friend __forceinline float reduce_add(const pairi& v) { return     reduce_add(v.l)+ reduce_add(v.h) ; }

    friend __forceinline size_t select_min(const pairi& v) { return __bsf(movemask(v == vreduce_min(v))); }
    friend __forceinline size_t select_max(const pairi& v) { return __bsf(movemask(v == vreduce_max(v))); }

    friend __forceinline size_t select_min(const Mask& valid, const pairi& v) { 
      const pairi a = select(valid,v,pairi(pos_inf)); 
      return __bsf(movemask(valid & (a == vreduce_min(a)))); 
    }
    friend __forceinline size_t select_max(const Mask& valid, const pairi& v) { 
      const pairi a = select(valid,v,pairi(neg_inf)); 
      return __bsf(movemask(valid & (a == vreduce_max(a)))); 
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Output Operators
    ////////////////////////////////////////////////////////////////////////////////
    
    friend inline std::ostream& operator<<(std::ostream& cout, const pairi& a) {
      return cout << "<" << a.l << ", " << a.h << ">";
    }
  };
}

#endif

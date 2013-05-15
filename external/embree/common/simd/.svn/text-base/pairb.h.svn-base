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

#ifndef __EMBREE_PAIRB_H__
#define __EMBREE_PAIRB_H__

namespace embree
{
  /*! Doubling SIMD bool types. */
  template<typename T> struct pairb
  {
    typedef pairb Mask;           // mask type for us
    enum   { size = 2*T::size };  // number of SIMD elements
    T l,h;                        // low and high data element
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Constructors, Assignment & Cast Operators
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline pairb           ( ) {}
    __forceinline pairb           ( const pairb& a ) { l = a.l; h = a.h; }
    __forceinline pairb& operator=( const pairb& a ) { l = a.l; h = a.h; return *this; }

    __forceinline pairb           ( const T& a             ) : l(a), h(a) {}
    __forceinline pairb           ( const T& a, const T& b ) : l(a), h(b) {}

    __forceinline pairb           ( bool  a ) : l(a), h(a) {}
    __forceinline pairb           ( bool  a, bool  b) : l(a), h(b) {}
    __forceinline pairb           ( bool  a, bool  b, bool  c, bool  d) : l(a,b), h(c,d) {}
    __forceinline pairb           ( bool  a, bool  b, bool  c, bool  d, bool  e, bool  f, bool  g, bool  h) : l(a,b,c,d), h(e,f,g,h) {}
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline pairb( FalseTy ) : l(False), h(False) {}
    __forceinline pairb( TrueTy  ) : l(True ), h(True ) {}

    ////////////////////////////////////////////////////////////////////////////////
    /// Array Access
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline const int32& operator []( const size_t i ) const { assert(i < size); return ((int32*)this)[i]; }
    __forceinline       int32& operator []( const size_t i )       { assert(i < size); return ((int32*)this)[i]; }

    ////////////////////////////////////////////////////////////////////////////////
    /// Unary Operators
    ////////////////////////////////////////////////////////////////////////////////
    
    friend __forceinline const pairb operator !( const pairb& a ) { return pairb(!a.l,!a.h); }
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Binary Operators
    ////////////////////////////////////////////////////////////////////////////////
    
    friend __forceinline const pairb operator &( const pairb& a, const pairb& b ) { return pairb(a.l & b.l, a.h & b.h); }
    friend __forceinline const pairb operator |( const pairb& a, const pairb& b ) { return pairb(a.l | b.l, a.h | b.h); }
    friend __forceinline const pairb operator ^( const pairb& a, const pairb& b ) { return pairb(a.l ^ b.l, a.h ^ b.h); }
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Assignment Operators
    ////////////////////////////////////////////////////////////////////////////////
    
    friend __forceinline pairb operator &=( pairb& a, const pairb& b ) { return a = a & b; }
    friend __forceinline pairb operator |=( pairb& a, const pairb& b ) { return a = a | b; }
    friend __forceinline pairb operator ^=( pairb& a, const pairb& b ) { return a = a ^ b; }
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Comparison Operators + Select
    ////////////////////////////////////////////////////////////////////////////////
    
    friend __forceinline const pairb operator !=( const pairb& a, const pairb& b ) { return pairb(a.l != b.l, a.h != b.h); }
    friend __forceinline const pairb operator ==( const pairb& a, const pairb& b ) { return pairb(a.l == b.l, a.h == b.h); }
    
    friend __forceinline const pairb select ( const Mask& m, const pairb& a, const pairb& b ) { 
      return pairb(select(m.l,a.l,b.l),select(m.h,a.h,b.h)); 
    }
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Reduction Operations
    ////////////////////////////////////////////////////////////////////////////////
    
    friend __forceinline size_t popcnt    ( const pairb& a ) { return popcnt    (a.l) + popcnt    (a.h); }
    friend __forceinline size_t movemask  ( const pairb& a ) { return (movemask(a.h) << T::size) | movemask(a.l); }
        
    friend __forceinline bool   reduce_and( const pairb& a ) { return reduce_and(a.l) & reduce_and(a.h); }
    friend __forceinline bool   reduce_or ( const pairb& a ) { return reduce_or (a.l) | reduce_or (a.h); }
    friend __forceinline bool   all       ( const pairb& a ) { return reduce_and(a); }
    friend __forceinline bool   any       ( const pairb& a ) { return reduce_or (a); }
    friend __forceinline bool   none      ( const pairb& a ) { return !reduce_or(a); }
    
    ////////////////////////////////////////////////////////////////////////////////
    /// Output Operators
    ////////////////////////////////////////////////////////////////////////////////
    
    friend inline std::ostream& operator<<(std::ostream& cout, const pairb& a) {
      return cout << "<" << a.l << ", " << a.h << ">";
    }
  };
}

#endif

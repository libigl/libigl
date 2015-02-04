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

#include "math/linearspace3.h"
#include "math/quaternion.h"
#include "math/bbox.h"

namespace embree
{
  #define VectorT typename L::Vector
  #define ScalarT typename L::Vector::Scalar

  ////////////////////////////////////////////////////////////////////////////////
  // Affine Space
  ////////////////////////////////////////////////////////////////////////////////

  template<typename L>
    struct AffineSpaceT
    {
      L l;           /*< linear part of affine space */
      VectorT p;     /*< affine part of affine space */

      ////////////////////////////////////////////////////////////////////////////////
      // Constructors, Assignment, Cast, Copy Operations
      ////////////////////////////////////////////////////////////////////////////////

      __forceinline AffineSpaceT           ( )                           { }
      __forceinline AffineSpaceT           ( const AffineSpaceT& other ) { l = other.l; p = other.p; }
      __forceinline AffineSpaceT           ( const L           & other ) { l = other  ; p = VectorT(zero); }
      __forceinline AffineSpaceT& operator=( const AffineSpaceT& other ) { l = other.l; p = other.p; return *this; }

      __forceinline AffineSpaceT( const VectorT& vx, const VectorT& vy, const VectorT& vz, const VectorT& p ) : l(vx,vy,vz), p(p) {}
      __forceinline AffineSpaceT( const L& l, const VectorT& p ) : l(l), p(p) {}

      template<typename L1> __forceinline explicit AffineSpaceT( const AffineSpaceT<L1>& s ) : l(s.l), p(s.p) {}

      ////////////////////////////////////////////////////////////////////////////////
      // Constants
      ////////////////////////////////////////////////////////////////////////////////

      __forceinline AffineSpaceT( ZeroTy ) : l(zero), p(zero) {}
      __forceinline AffineSpaceT( OneTy )  : l(one),  p(zero) {}

      /*! return matrix for scaling */
      static __forceinline AffineSpaceT scale(const VectorT& s) { return AffineSpaceT(L::scale(s),zero); }

      /*! return matrix for translation */
      static __forceinline AffineSpaceT translate(const VectorT& p) { return AffineSpaceT(one,p); }

      /*! return matrix for rotation around arbitrary axis */
      static __forceinline AffineSpaceT rotate(const VectorT& u, const ScalarT& r) { return AffineSpaceT(L::rotate(u,r),zero); }

      /*! return matrix for rotation around arbitrary axis and point */
      static __forceinline AffineSpaceT rotate(const VectorT& p, const VectorT& u, const ScalarT& r) { return translate(+p) * rotate(u,r) * translate(-p);  }

      /*! return matrix for looking at given point */
      static __forceinline AffineSpaceT lookat(const VectorT& eye, const VectorT& point, const VectorT& up) {
        VectorT Z = normalize(point-eye);
        VectorT U = normalize(cross(up,Z));
        VectorT V = normalize(cross(Z,U));
        return AffineSpaceT(L(U,V,Z),eye);
      }

    };

  ////////////////////////////////////////////////////////////////////////////////
  // Unary Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename L> __forceinline AffineSpaceT<L> operator -( const AffineSpaceT<L>& a ) { return AffineSpaceT<L>(-a.l,-a.p); }
  template<typename L> __forceinline AffineSpaceT<L> operator +( const AffineSpaceT<L>& a ) { return AffineSpaceT<L>(+a.l,+a.p); }
  template<typename L> __forceinline AffineSpaceT<L>        rcp( const AffineSpaceT<L>& a ) { L il = rcp(a.l); return AffineSpaceT<L>(il,-(il*a.p)); }

  ////////////////////////////////////////////////////////////////////////////////
  // Binary Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename L> __forceinline const AffineSpaceT<L> operator +( const AffineSpaceT<L>& a, const AffineSpaceT<L>& b ) { return AffineSpaceT<L>(a.l+b.l,a.p+b.p); }
  template<typename L> __forceinline const AffineSpaceT<L> operator -( const AffineSpaceT<L>& a, const AffineSpaceT<L>& b ) { return AffineSpaceT<L>(a.l-b.l,a.p-b.p); }

  template<typename L> __forceinline const AffineSpaceT<L> operator *( const ScalarT        & a, const AffineSpaceT<L>& b ) { return AffineSpaceT<L>(a*b.l,a*b.p); }
  template<typename L> __forceinline const AffineSpaceT<L> operator *( const AffineSpaceT<L>& a, const AffineSpaceT<L>& b ) { return AffineSpaceT<L>(a.l*b.l,a.l*b.p+a.p); }
  template<typename L> __forceinline const AffineSpaceT<L> operator /( const AffineSpaceT<L>& a, const AffineSpaceT<L>& b ) { return a * rcp(b); }
  template<typename L> __forceinline const AffineSpaceT<L> operator /( const AffineSpaceT<L>& a, const ScalarT        & b ) { return a * rcp(b); }

  template<typename L> __forceinline AffineSpaceT<L>& operator *=( AffineSpaceT<L>& a, const AffineSpaceT<L>& b ) { return a = a * b; }
  template<typename L> __forceinline AffineSpaceT<L>& operator *=( AffineSpaceT<L>& a, const ScalarT        & b ) { return a = a * b; }
  template<typename L> __forceinline AffineSpaceT<L>& operator /=( AffineSpaceT<L>& a, const AffineSpaceT<L>& b ) { return a = a / b; }
  template<typename L> __forceinline AffineSpaceT<L>& operator /=( AffineSpaceT<L>& a, const ScalarT        & b ) { return a = a / b; }

  template<typename L> __forceinline const VectorT xfmPoint (const AffineSpaceT<L>& m, const VectorT& p) { return madd(VectorT(p.x),m.l.vx,madd(VectorT(p.y),m.l.vy,madd(VectorT(p.z),m.l.vz,m.p))); }
  template<typename L> __forceinline const VectorT xfmVector(const AffineSpaceT<L>& m, const VectorT& v) { return xfmVector(m.l,v); }
  template<typename L> __forceinline const VectorT xfmNormal(const AffineSpaceT<L>& m, const VectorT& n) { return xfmNormal(m.l,n); }

  __forceinline const BBox<Vec3fa> xfmBounds(const AffineSpaceT<LinearSpace3<Vec3fa> >& m, const BBox<Vec3fa>& b) 
  { 
    BBox3fa dst = empty;
    const Vec3fa p0(b.lower.x,b.lower.y,b.lower.z); dst.extend(xfmPoint(m,p0));
    const Vec3fa p1(b.lower.x,b.lower.y,b.upper.z); dst.extend(xfmPoint(m,p1));
    const Vec3fa p2(b.lower.x,b.upper.y,b.lower.z); dst.extend(xfmPoint(m,p2));
    const Vec3fa p3(b.lower.x,b.upper.y,b.upper.z); dst.extend(xfmPoint(m,p3));
    const Vec3fa p4(b.upper.x,b.lower.y,b.lower.z); dst.extend(xfmPoint(m,p4));
    const Vec3fa p5(b.upper.x,b.lower.y,b.upper.z); dst.extend(xfmPoint(m,p5));
    const Vec3fa p6(b.upper.x,b.upper.y,b.lower.z); dst.extend(xfmPoint(m,p6));
    const Vec3fa p7(b.upper.x,b.upper.y,b.upper.z); dst.extend(xfmPoint(m,p7));
    return dst;
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Comparison Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename L> __forceinline bool operator ==( const AffineSpaceT<L>& a, const AffineSpaceT<L>& b ) { return a.l == b.l && a.p == b.p; }
  template<typename L> __forceinline bool operator !=( const AffineSpaceT<L>& a, const AffineSpaceT<L>& b ) { return a.l != b.l || a.p != b.p; }

  ////////////////////////////////////////////////////////////////////////////////
  // Output Operators
  ////////////////////////////////////////////////////////////////////////////////

  template<typename L> static std::ostream& operator<<(std::ostream& cout, const AffineSpaceT<L>& m) {
    return cout << "{ l = " << m.l << ", p = " << m.p << " }";
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Template Instantiations
  ////////////////////////////////////////////////////////////////////////////////

  typedef AffineSpaceT<LinearSpace3f> AffineSpace3f;
  typedef AffineSpaceT<LinearSpace3fa> AffineSpace3fa;
  typedef AffineSpaceT<Quaternion3f > OrthonormalSpace3f;

  #undef VectorT
  #undef ScalarT
}

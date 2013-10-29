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

#ifndef __EMBREE_AFFINE_SPACE_H__
#define __EMBREE_AFFINE_SPACE_H__

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
      static __forceinline AffineSpaceT lookAtPoint(const VectorT& eye, const VectorT& point, const VectorT& up) {
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

  template<typename L> __forceinline const VectorT xfmPoint (const AffineSpaceT<L>& m, const VectorT& p) { return xfmPoint(m.l,p) + m.p; }
  template<typename L> __forceinline const VectorT xfmVector(const AffineSpaceT<L>& m, const VectorT& v) { return xfmVector(m.l,v); }
  template<typename L> __forceinline const VectorT xfmNormal(const AffineSpaceT<L>& m, const VectorT& n) { return xfmNormal(m.l,n); }

  __forceinline const BBox<Vector3f> xfmBounds(const AffineSpaceT<LinearSpace3<Vector3f> >& m, const BBox<Vector3f>& b) 
  { 
    BBox3f dst = empty;
    const Vector3f p0(b.lower.x,b.lower.y,b.lower.z); dst.grow(xfmPoint(m,p0));
    const Vector3f p1(b.lower.x,b.lower.y,b.upper.z); dst.grow(xfmPoint(m,p1));
    const Vector3f p2(b.lower.x,b.upper.y,b.lower.z); dst.grow(xfmPoint(m,p2));
    const Vector3f p3(b.lower.x,b.upper.y,b.upper.z); dst.grow(xfmPoint(m,p3));
    const Vector3f p4(b.upper.x,b.lower.y,b.lower.z); dst.grow(xfmPoint(m,p4));
    const Vector3f p5(b.upper.x,b.lower.y,b.upper.z); dst.grow(xfmPoint(m,p5));
    const Vector3f p6(b.upper.x,b.upper.y,b.lower.z); dst.grow(xfmPoint(m,p6));
    const Vector3f p7(b.upper.x,b.upper.y,b.upper.z); dst.grow(xfmPoint(m,p7));
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
  //typedef AffineSpaceT<LinearSpace3d> AffineSpace3d;
  typedef AffineSpaceT<Quaternion3f > OrthonormalSpace3f;

  ////////////////////////////////////////////////////////////////////////////////
  // Data conversions for AffineSpace3f
  ////////////////////////////////////////////////////////////////////////////////

  struct Array12f {
    float values[12];
    operator float*() { return(values); }
  };
  
  __forceinline Array12f copyToArray(const AffineSpace3f& xfm)  
  {
    Array12f values;
    values[ 0] = xfm.l.vx.x;  values[ 1] = xfm.l.vx.y;  values[ 2] = xfm.l.vx.z;       
    values[ 3] = xfm.l.vy.x;  values[ 4] = xfm.l.vy.y;  values[ 5] = xfm.l.vy.z;       
    values[ 6] = xfm.l.vz.x;  values[ 7] = xfm.l.vz.y;  values[ 8] = xfm.l.vz.z;       
    values[ 9] = xfm.p.x;     values[10] = xfm.p.y;     values[11] = xfm.p.z;       
    return values;
  }
  
  __forceinline AffineSpace3f copyFromArray(const float* v) 
  {
    return AffineSpace3f(LinearSpace3f(Vector3f(v[0],v[1],v[2]),
                                       Vector3f(v[3],v[4],v[5]),
                                       Vector3f(v[6],v[7],v[8])),
                         Vector3f(v[9],v[10],v[11]));
  }

  #undef VectorT
  #undef ScalarT
}

#endif

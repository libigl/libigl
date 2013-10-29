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

#ifndef __EMBREE_BBOX_H__
#define __EMBREE_BBOX_H__

#include "vec2.h"
#include "vec3.h"

namespace embree
{
  template<typename T>
  struct BBox
  {
    T lower, upper;

    ////////////////////////////////////////////////////////////////////////////////
    /// Construction
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline BBox           ( )                   { }
    __forceinline BBox           ( const BBox& other ) { lower = other.lower; upper = other.upper; }
    __forceinline BBox& operator=( const BBox& other ) { lower = other.lower; upper = other.upper; return *this; }

    __forceinline BBox ( const T& v                 ) : lower(v),   upper(v) {}
    __forceinline BBox ( const T& lower, const T& upper ) : lower(lower), upper(upper) {}

    ////////////////////////////////////////////////////////////////////////////////
    /// Growing Bounds
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline void grow(const BBox& other) { lower = min(lower,other.lower); upper = max(upper,other.upper); }
    __forceinline void grow(const T   & other) { lower = min(lower,other    ); upper = max(upper,other    ); }

    /*! computes the size of the box */
    __forceinline bool empty() const { for (int i=0; i<T::N; i++) if (lower[i] > upper[i]) return true; return false; }

    /*! computes the size of the box */
    __forceinline T size() const { return upper - lower; }

    ////////////////////////////////////////////////////////////////////////////////
    /// Constants
    ////////////////////////////////////////////////////////////////////////////////

    __forceinline BBox( EmptyTy ) : lower(pos_inf), upper(neg_inf) {}
    __forceinline BBox( FullTy  ) : lower(neg_inf), upper(pos_inf) {}
    __forceinline BBox( FalseTy ) : lower(pos_inf), upper(neg_inf) {}
    __forceinline BBox( TrueTy  ) : lower(neg_inf), upper(pos_inf) {}
    __forceinline BBox( NegInfTy ): lower(pos_inf), upper(neg_inf) {}
    __forceinline BBox( PosInfTy ): lower(neg_inf), upper(pos_inf) {}
  };

  /*! tests if box is empty */
  template<typename T> __forceinline bool isEmpty(const BBox<T>& box) { for (int i=0; i<T::N; i++) if (box.lower[i] > box.upper[i]) return true; return false; }

  /*! computes the center of the box */
  template<typename T> __forceinline const T center (const BBox<T>& box) { return T(.5f)*(box.lower + box.upper); }
  template<typename T> __forceinline const T center2(const BBox<T>& box) { return box.lower + box.upper; }

  /*! computes the size of the box */
  template<typename T> __forceinline const T size(const BBox<T>& box) { return box.upper - box.lower; }

  /*! computes the volume of a bounding box */
  template<typename T> __forceinline const T volume( const BBox<T>& b ) { return reduce_mul(size(b)); }

  /*! computes the surface area of a bounding box */
  template<typename T> __forceinline const T     area( const BBox<Vec2<T> >& b ) { const Vec2<T> d = size(b); return d.x*d.y; }
  __forceinline float     area( const BBox<Vector3f>& b ) { const Vector3f d = size(b); return 2.0f*(d.x*(d.y+d.z)+d.y*d.z); }
  __forceinline float halfArea( const BBox<Vector3f>& b ) { const Vector3f d = size(b); return d.x*(d.y+d.z)+d.y*d.z; }

  /*! merges bounding boxes and points */
  template<typename T> __forceinline const BBox<T> merge( const BBox<T>& a, const       T& b ) { return BBox<T>(min(a.lower, b    ), max(a.upper, b    )); }
  template<typename T> __forceinline const BBox<T> merge( const       T& a, const BBox<T>& b ) { return BBox<T>(min(a    , b.lower), max(a    , b.upper)); }
  template<typename T> __forceinline const BBox<T> merge( const BBox<T>& a, const BBox<T>& b ) { return BBox<T>(min(a.lower, b.lower), max(a.upper, b.upper)); }
  template<typename T> __forceinline const BBox<T> merge( const BBox<T>& a, const BBox<T>& b, const BBox<T>& c ) { return merge(a,merge(b,c)); }
  template<typename T> __forceinline const BBox<T>& operator+=( BBox<T>& a, const BBox<T>& b ) { return a = merge(a,b); }
  template<typename T> __forceinline const BBox<T>& operator+=( BBox<T>& a, const       T& b ) { return a = merge(a,b); }

  /*! Merges four boxes. */
  template<typename T> __forceinline BBox<T> merge(const BBox<T>& a, const BBox<T>& b, const BBox<T>& c, const BBox<T>& d) {
    return merge(merge(a,b),merge(c,d));
  }

  /*! Merges eight boxes. */
  template<typename T> __forceinline BBox<T> merge(const BBox<T>& a, const BBox<T>& b, const BBox<T>& c, const BBox<T>& d,
                                                   const BBox<T>& e, const BBox<T>& f, const BBox<T>& g, const BBox<T>& h) {
    return merge(merge(a,b,c,d),merge(e,f,g,h));
  }

  /*! Comparison Operators */
  template<typename T> __forceinline bool operator==( const BBox<T>& a, const BBox<T>& b ) { return a.lower == b.lower && a.upper == b.upper; }
  template<typename T> __forceinline bool operator!=( const BBox<T>& a, const BBox<T>& b ) { return a.lower != b.lower || a.upper != b.upper; }

  /*! scaling */
  template<typename T> __forceinline BBox<T> operator *( const float& a, const BBox<T>& b ) { return BBox<T>(a*b.lower,a*b.upper); }

  /*! intersect bounding boxes */
  template<typename T> __forceinline const BBox<T> intersect( const BBox<T>& a, const BBox<T>& b ) { return BBox<T>(max(a.lower, b.lower), min(a.upper, b.upper)); }
  template<typename T> __forceinline const BBox<T> intersect( const BBox<T>& a, const BBox<T>& b, const BBox<T>& c ) { return intersect(a,intersect(b,c)); }

  /*! tests if bounding boxes (and points) are disjoint (empty intersection) */
  template<typename T> __inline bool disjoint( const BBox<T>& a, const BBox<T>& b )
  { const T d = min(a.upper, b.upper) - max(a.lower, b.lower); for ( size_t i = 0 ; i < T::N ; i++ ) if ( d[i] < typename T::Scalar(zero) ) return true; return false; }
  template<typename T> __inline bool disjoint( const BBox<T>& a, const  T& b )
  { const T d = min(a.upper, b)       - max(a.lower, b);       for ( size_t i = 0 ; i < T::N ; i++ ) if ( d[i] < typename T::Scalar(zero) ) return true; return false; }
  template<typename T> __inline bool disjoint( const  T& a, const BBox<T>& b )
  { const T d = min(a, b.upper)       - max(a, b.lower);       for ( size_t i = 0 ; i < T::N ; i++ ) if ( d[i] < typename T::Scalar(zero) ) return true; return false; }

  /*! tests if bounding boxes (and points) are conjoint (non-empty intersection) */
  template<typename T> __inline bool conjoint( const BBox<T>& a, const BBox<T>& b )
  { const T d = min(a.upper, b.upper) - max(a.lower, b.lower); for ( size_t i = 0 ; i < T::N ; i++ ) if ( d[i] < typename T::Scalar(zero) ) return false; return true; }
  template<typename T> __inline bool conjoint( const BBox<T>& a, const  T& b )
  { const T d = min(a.upper, b)       - max(a.lower, b);       for ( size_t i = 0 ; i < T::N ; i++ ) if ( d[i] < typename T::Scalar(zero) ) return false; return true; }
  template<typename T> __inline bool conjoint( const  T& a, const BBox<T>& b )
  { const T d = min(a, b.upper)       - max(a, b.lower);       for ( size_t i = 0 ; i < T::N ; i++ ) if ( d[i] < typename T::Scalar(zero) ) return false; return true; }

  /*! subset relation */
  template<typename T> __inline bool subset( const BBox<T>& a, const BBox<T>& b )
  { 
    for ( size_t i = 0 ; i < T::N ; i++ ) if ( a.lower[i] < b.lower[i] ) return false;
    for ( size_t i = 0 ; i < T::N ; i++ ) if ( a.upper[i] > b.upper[i] ) return false;
    return true; 
  }

  /*! output operator */
  template<typename T> __forceinline std::ostream& operator<<(std::ostream& cout, const BBox<T>& box) {
    return cout << "[" << box.lower << "; " << box.upper << "]";
  }

  /*! default template instantiations */
  typedef BBox<Vec2f > BBox2f;
  typedef BBox<Vector3f > BBox3f;
}

#endif

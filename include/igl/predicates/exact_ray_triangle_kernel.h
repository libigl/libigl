// This file is part of libigl, a simple c++ geometry processing library.
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PREDICATES_EXACT_RAY_TRIANGLE_KERNEL_H
#define IGL_PREDICATES_EXACT_RAY_TRIANGLE_KERNEL_H

#include "exact_point_triangle_kernel.h"

namespace igl { namespace predicates { namespace exact_ray_triangle_kernel
{
  namespace k=exact_point_triangle_kernel;
  using Point=k::Point;
  template<typename Number> inline auto direction_vector(const Point &direction)
  {
    return k::Vector<k::Input<Number>>{{direction[0],direction[1],direction[2]}};
  }
  template<typename Number> inline auto normal(const Point &a,const Point &b,const Point &c)
  { return k::cross(k::difference<Number>(b,a),k::difference<Number>(c,a)); }

  // -1 requests a more accurate stage; 0/1/2 are None/Hit/Coplanar.
  template<typename Number> inline int classify(const Point &origin,const Point &direction,
    const Point &a,const Point &b,const Point &c)
  {
    const auto ray=direction_vector<Number>(direction);
    const auto n=normal<Number>(a,b,c);
    const int denominator=k::sign(k::dot(ray,n));
    if(denominator==2) return -1;
    const auto ao=k::difference<Number>(a,origin);
    const int numerator=[&]()
    {
      if constexpr(std::is_same<Number,k::Expansion>::value)
        return k::orientation(a,b,c,origin);
      else return k::sign(k::dot(n,ao));
    }();
    if(numerator==2) return -1;
    if(denominator==0) return numerator==0 ? 2 : 0;
    if(numerator==-denominator) return 0;
    const auto bo=k::difference<Number>(b,origin);
    const int first=k::sign(k::dot(ray,k::cross(ao,bo)));
    if(first==2) return -1;
    if(first==-denominator) return 0;
    const auto co=k::difference<Number>(c,origin);
    const int second=k::sign(k::dot(ray,k::cross(bo,co)));
    if(second==2) return -1;
    if(second==-denominator) return 0;
    const int third=k::sign(k::dot(ray,k::cross(co,ao)));
    if(third==2) return -1;
    return third==-denominator ? 0 : 1;
  }
  inline int classify(const Point &origin,const Point &direction,
    const Point &a,const Point &b,const Point &c)
  {
    const int fast=classify<k::Filter>(origin,direction,a,b,c);
    if(fast>=0) return fast;
    // Exact classification is cheap enough that retrying every uncertain
    // determinant with nextafter intervals costs more on thin triangles.
    std::array<Point,4> scaled{{origin,a,b,c}};
    std::array<Point,1> ray{{direction}};
    k::rescale<3>(scaled); k::rescale<3>(ray);
    return classify<k::Expansion>(scaled[0],ray[0],scaled[1],scaled[2],scaled[3]);
  }

  template<typename Number> inline int compare(const Point &origin,const Point &direction,
    const Point &a1,const Point &b1,const Point &c1,
    const Point &a2,const Point &b2,const Point &c2)
  {
    const auto ray=direction_vector<Number>(direction);
    const auto normal1=normal<Number>(a1,b1,c1),normal2=normal<Number>(a2,b2,c2);
    const auto denominator1=k::dot(ray,normal1),denominator2=k::dot(ray,normal2);
    const int sign1=k::sign(denominator1),sign2=k::sign(denominator2);
    if(sign1==2 || sign2==2) return 2;
    if(sign1==0 || sign2==0)
      throw std::invalid_argument("Hit ordering requires two transverse triangles");
    const auto numerator1=k::dot(normal1,k::difference<Number>(a1,origin));
    const auto numerator2=k::dot(normal2,k::difference<Number>(a2,origin));
    const int order=k::sign(numerator1*denominator2-numerator2*denominator1);
    // Dividing by denominators with opposite signs reverses the comparison.
    return order==2 ? 2 : order*sign1*sign2;
  }
  inline int compare(const Point &origin,const Point &direction,
    const Point &a1,const Point &b1,const Point &c1,
    const Point &a2,const Point &b2,const Point &c2)
  {
    const int fast=compare<k::Filter>(origin,direction,a1,b1,c1,a2,b2,c2);
    if(fast!=2) return fast;
    const int interval=compare<k::Interval>(origin,direction,a1,b1,c1,a2,b2,c2);
    if(interval!=2) return interval;
    std::array<Point,7> scaled{{origin,a1,b1,c1,a2,b2,c2}};
    std::array<Point,1> ray{{direction}};
    k::rescale<6>(scaled); k::rescale<6>(ray);
    return compare<k::Expansion>(scaled[0],ray[0],scaled[1],scaled[2],scaled[3],
      scaled[4],scaled[5],scaled[6]);
  }
}}}
#endif

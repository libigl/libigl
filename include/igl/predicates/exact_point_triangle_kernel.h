// This file is part of libigl, a simple c++ geometry processing library.
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PREDICATES_EXACT_POINT_TRIANGLE_KERNEL_H
#define IGL_PREDICATES_EXACT_POINT_TRIANGLE_KERNEL_H
#include <predicates.h>
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace igl { namespace predicates { namespace exact_point_triangle_kernel
{
  static_assert(std::is_same<IGL_PREDICATES_REAL,double>::value,
    "Exact metric predicates require the double-precision predicates library");
  // Private arithmetic objects never escape a predicate invocation.
  constexpr int capacity = std::numeric_limits<double>::max_exponent
    - std::numeric_limits<double>::min_exponent + std::numeric_limits<double>::digits;
  struct Expansion
  {
    int size;
    double data[capacity];
    Expansion(double value=0):size(1) { data[0]=value; }
    Expansion(const Expansion &other):size(other.size)
    { std::copy_n(other.data,size,data); }
    Expansion &operator=(const Expansion &other)
    { size=other.size; std::copy_n(other.data,size,data); return *this; }
    void assign(int length,const double *values)
    {
      if(length>capacity || !std::isfinite(values[length-1]))
        throw std::overflow_error("Exact metric expansion exceeded its arithmetic range");
      size=length; std::copy_n(values,size,data);
    }
  };
  inline int sign(const Expansion &value)
  { return (value.data[value.size-1]>0)-(value.data[value.size-1]<0); }
  inline Expansion operator-(const Expansion &a,const Expansion &b)
  {
    Expansion result;
    if(a.size==1 && b.size==1)
    {
      double head,tail;
      ::predicates_two_diff(a.data[0],b.data[0],&head,&tail);
      result.data[0]=tail==0 ? head : tail;
      if(tail!=0 && head!=0) { result.data[1]=head; result.size=2; }
      return result;
    }
    double negative[capacity],sum[2*capacity];
    for(int i=0;i<b.size;++i) negative[i]=-b.data[i];
    result.assign(::predicates_expansion_sum(a.size,a.data,b.size,negative,sum),sum);
    return result;
  }
  inline Expansion operator+(const Expansion &a,const Expansion &b)
  {
    Expansion result;
    double sum[2*capacity];
    result.assign(::predicates_expansion_sum(a.size,a.data,b.size,b.data,sum),sum);
    return result;
  }
  inline Expansion operator*(const Expansion &a,const Expansion &b)
  {
    Expansion result;
    if(sign(a)==0 || sign(b)==0) return result;
    double scaled[2*capacity],sum[3*capacity];
    for(int i=0;i<b.size;++i)
    {
      const int length=::predicates_expansion_scale(a.size,a.data,b.data[i],scaled);
      result.assign(::predicates_expansion_sum(result.size,result.data,length,scaled,sum),sum);
    }
    return result;
  }

  // The largest expression below has degree ten. Counting rounding factors
  // in each expanded monomial gives at most 30 (including coordinate
  // differences), both for its value and for this absolute permanent.
  // 256*epsilon safely exceeds gamma_30/(1-gamma_30) in round-to-nearest.
  // See experiments/exact-points/README.md for the recurrence and assumptions.
  struct Filter
  {
    double value,permanent;
    Filter(double x=0):value(x),permanent(std::abs(x)) {}
    Filter(double x,double bound):value(x),permanent(bound) {}
  };
  inline Filter operator+(const Filter &a,const Filter &b)
  { return {a.value+b.value,a.permanent+b.permanent}; }
  inline Filter operator-(const Filter &a,const Filter &b)
  { return {a.value-b.value,a.permanent+b.permanent}; }
  inline Filter operator*(const Filter &a,const Filter &b)
  { return {a.value*b.value,a.permanent*b.permanent}; }
  inline int sign(const Filter &a)
  {
    const double error=256*std::numeric_limits<double>::epsilon()*a.permanent;
    if(a.value>error) return 1;
    if(a.value<-error) return -1;
    return a.value==0 && a.permanent==0 ? 0 : 2;
  }

  // Outward rounding without changing the thread's rounding mode. Unknown
  // signs (2) trigger a fresh evaluation using expansions, never a guess.
  struct Interval
  {
    double lo,hi;
    Interval(double value=0):lo(value),hi(value) {}
    Interval(double lower,double upper):lo(lower),hi(upper) {}
  };
  inline double down(double x)
  { return std::nextafter(x,-std::numeric_limits<double>::infinity()); }
  inline double up(double x)
  { return std::nextafter(x,std::numeric_limits<double>::infinity()); }
  inline int sign(const Interval &value)
  {
    if(value.lo>0) return 1;
    if(value.hi<0) return -1;
    return value.lo==0 && value.hi==0 ? 0 : 2;
  }
  inline Interval operator+(const Interval &a,const Interval &b)
  {
    if(a.lo==0 && a.hi==0) return b;
    if(b.lo==0 && b.hi==0) return a;
    return {down(a.lo+b.lo),up(a.hi+b.hi)};
  }
  inline Interval operator-(const Interval &a,const Interval &b)
  {
    if(a.lo==a.hi && b.lo==b.hi && a.lo==b.lo) return Interval(0);
    if(b.lo==0 && b.hi==0) return a;
    return {down(a.lo-b.hi),up(a.hi-b.lo)};
  }
  inline Interval operator*(const Interval &a,const Interval &b)
  {
    if((a.lo==0 && a.hi==0)||(b.lo==0 && b.hi==0)) return Interval(0);
    const double products[4]={a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};
    for(double value:products) if(std::isnan(value))
      return {-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()};
    return {down(*std::min_element(products,products+4)),
            up(*std::max_element(products,products+4))};
  }
  using Point=std::array<double,3>;
  template<typename Derived> inline Point point(const Eigen::MatrixBase<Derived> &p)
  {
    static_assert(std::is_same<typename Derived::Scalar,double>::value ||
      std::is_same<typename Derived::Scalar,float>::value,"Exact metric inputs must be float or double");
    if(p.size()!=3) throw std::invalid_argument("Expected a three-coordinate point");
    Point result;
    for(int i=0;i<3;++i)
    {
      result[i]=p.derived().coeff(i);
      if(!std::isfinite(result[i])) throw std::invalid_argument("Exact metric inputs must be finite");
    }
    return result;
  }
  template<typename Number> using Vector=std::array<Number,3>;
  template<typename Number> inline Vector<Number> difference(const Point &a,const Point &b)
  { return {{Number(a[0])-Number(b[0]),Number(a[1])-Number(b[1]),Number(a[2])-Number(b[2])}}; }
  template<> inline Vector<Filter> difference<Filter>(const Point &a,const Point &b)
  { return {{Filter(a[0]-b[0]),Filter(a[1]-b[1]),Filter(a[2]-b[2])}}; }
  template<typename Number> inline Number dot(const Vector<Number> &a,const Vector<Number> &b)
  { return (a[0]*b[0]+a[1]*b[1])+a[2]*b[2]; }
  template<typename Number> inline Vector<Number> cross(const Vector<Number> &a,const Vector<Number> &b)
  { return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}}; }

  // Dimension 0/1/2 and corner/starting-corner/0. -1 means uncertain.
  template<typename Number> inline int feature(const Point &q,
    const std::array<Point,3> &p,int &id)
  {
    for(int i=0;i<3;++i)
    {
      const auto offset=difference<Number>(q,p[i]);
      const int a=sign(dot(offset,difference<Number>(p[(i+1)%3],p[i])));
      const int b=sign(dot(offset,difference<Number>(p[(i+2)%3],p[i])));
      if(a==1 || b==1) continue;
      if(a==2 || b==2) return -1;
      id=i; return 0;
    }
    for(int i=0;i<3;++i)
    {
      const auto edge=difference<Number>(p[(i+1)%3],p[i]);
      const auto inward=difference<Number>(p[(i+2)%3],p[i]);
      const auto offset=difference<Number>(q,p[i]);
      const auto projection=dot(offset,edge),length=dot(edge,edge);
      const int first=sign(projection),last=sign(length-projection);
      if(first==2 || last==2) return -1;
      if(first<=0 || last<=0) continue;
      const int side=sign(length*dot(offset,inward)-projection*dot(edge,inward));
      if(side==2) return -1;
      if(side<=0) { id=i; return 1; }
    }
    id=0; return 2;
  }
  inline int feature(const Point &q,const std::array<Point,3> &p,int &id)
  {
    const int fast=feature<Filter>(q,p,id);
    if(fast>=0) return fast;
    const int filtered=feature<Interval>(q,p,id);
    return filtered<0 ? feature<Expansion>(q,p,id) : filtered;
  }
  template<typename Number> struct Fraction { Number numerator,denominator; };
  template<typename Number> inline Fraction<Number> distance(const Point &q,
    const std::array<Point,3> &p,int dimension,int id)
  {
    const auto offset=difference<Number>(q,p[id]);
    if(dimension==0) return {dot(offset,offset),Number(1)};
    const auto edge=difference<Number>(p[(id+1)%3],p[id]);
    if(dimension==1)
    {
      const auto length=dot(edge,edge),projection=dot(offset,edge);
      return {dot(offset,offset)*length-projection*projection,length};
    }
    const auto normal=cross(edge,difference<Number>(p[(id+2)%3],p[id]));
    const auto height=dot(offset,normal);
    return {height*height,dot(normal,normal)};
  }
  template<typename Number> inline int compare(const Point &q,
    const std::array<Point,3> &a,int ad,int ai,
    const std::array<Point,3> &b,int bd,int bi)
  {
    if(ad==0 && bd==0 && a[ai]==b[bi]) return 0;
    if(ad==1 && bd==1 &&
       ((a[ai]==b[bi] && a[(ai+1)%3]==b[(bi+1)%3]) ||
        (a[ai]==b[(bi+1)%3] && a[(ai+1)%3]==b[bi]))) return 0;
    const auto first=distance<Number>(q,a,ad,ai),second=distance<Number>(q,b,bd,bi);
    return sign(first.numerator*second.denominator-second.numerator*first.denominator);
  }
  // Reuse exact feature classifications within a query. No constructed point
  // or expansion survives a predicate call; only dimension and local id do.
  inline int compare(const Point &q,const std::array<Point,3> &a,int ad,int ai,
    const std::array<Point,3> &b,int bd,int bi)
  {
    const int fast=compare<Filter>(q,a,ad,ai,b,bd,bi);
    if(fast!=2) return fast;
    const int filtered=compare<Interval>(q,a,ad,ai,b,bd,bi);
    return filtered==2 ? compare<Expansion>(q,a,ad,ai,b,bd,bi) : filtered;
  }
  template<typename Number> inline int compare(const Point &q,
    const std::array<Point,3> &a,const std::array<Point,3> &b)
  {
    int ai,bi;
    const int ad=feature<Number>(q,a,ai),bd=feature<Number>(q,b,bi);
    if(ad<0 || bd<0) return 2;
    return compare<Number>(q,a,ad,ai,b,bd,bi);
  }
  inline int compare(const Point &q,const std::array<Point,3> &a,const std::array<Point,3> &b)
  {
    const int fast=compare<Filter>(q,a,b);
    if(fast!=2) return fast;
    const int filtered=compare<Interval>(q,a,b);
    return filtered==2 ? compare<Expansion>(q,a,b) : filtered;
  }
  template<typename Number> inline int perpendicular_compare(const Point &q,
    const Point &v,const Point &a,const Point &b)
  {
    const auto offset=difference<Number>(q,v),first=difference<Number>(a,v),second=difference<Number>(b,v);
    const auto x=dot(offset,first),y=dot(offset,second);
    return sign(x*x*dot(second,second)-y*y*dot(first,first));
  }
  inline int perpendicular_compare(const Point &q,const Point &v,const Point &a,const Point &b)
  {
    const int fast=perpendicular_compare<Filter>(q,v,a,b);
    if(fast!=2) return fast;
    const int filtered=perpendicular_compare<Interval>(q,v,a,b);
    return filtered==2 ? perpendicular_compare<Expansion>(q,v,a,b) : filtered;
  }
  inline int orientation(const Point &a,const Point &b,const Point &c,const Point &d)
  {
    // The C API predates const-correctness; it does not modify these copies.
    Point pa=a,pb=b,pc=c,pd=d;
    const double value=::orient3d(pa.data(),pb.data(),pc.data(),pd.data());
    return (value>0)-(value<0);
  }
}}}
#endif

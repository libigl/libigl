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
  constexpr int bounded_capacity(int size) { return size<capacity ? size : capacity; }
  // Track expression-length bounds so coordinate differences and normals do
  // not each reserve the maximum binary64 expansion. Arithmetic scratch still
  // accommodates the full sum/scale output before checking its finite bound.
  template<int Capacity> struct ExpansionStorage
  {
    static_assert(Capacity>0 && Capacity<=capacity,"Invalid expansion capacity");
    int size;
    double data[Capacity];
    ExpansionStorage(double value=0):size(1) { data[0]=value; }
    ExpansionStorage(const ExpansionStorage &other):size(other.size)
    { std::copy_n(other.data,size,data); }
    ExpansionStorage &operator=(const ExpansionStorage &other)
    { size=other.size; std::copy_n(other.data,size,data); return *this; }
    template<int Other> ExpansionStorage(const ExpansionStorage<Other> &other):size(other.size)
    {
      static_assert(Other<=Capacity,"An expansion conversion must preserve its capacity");
      std::copy_n(other.data,size,data);
    }
    template<int Other> ExpansionStorage &operator=(const ExpansionStorage<Other> &other)
    {
      static_assert(Other<=Capacity,"An expansion assignment must preserve its capacity");
      size=other.size; std::copy_n(other.data,size,data); return *this;
    }
    void assign(int length,const double *values)
    {
      if(length>Capacity || !std::isfinite(values[length-1]))
        throw std::overflow_error("Exact metric expansion exceeded its arithmetic range");
      size=length; std::copy_n(values,size,data);
    }
  };
  using Expansion=ExpansionStorage<capacity>;
  template<int Capacity> inline int sign(const ExpansionStorage<Capacity> &value)
  { return (value.data[value.size-1]>0)-(value.data[value.size-1]<0); }
  template<int A,int B> inline ExpansionStorage<bounded_capacity(A+B)> operator-(
    const ExpansionStorage<A> &a,const ExpansionStorage<B> &b)
  {
    ExpansionStorage<bounded_capacity(A+B)> result;
    if(a.size==1 && b.size==1)
    {
      double head,tail;
      ::predicates_two_diff(a.data[0],b.data[0],&head,&tail);
      result.data[0]=tail==0 ? head : tail;
      if(tail!=0 && head!=0) { result.data[1]=head; result.size=2; }
      return result;
    }
    double negative[B],sum[A+B];
    for(int i=0;i<b.size;++i) negative[i]=-b.data[i];
    result.assign(::predicates_expansion_sum(a.size,a.data,b.size,negative,sum),sum);
    return result;
  }
  template<int A,int B> inline ExpansionStorage<bounded_capacity(A+B)> operator+(
    const ExpansionStorage<A> &a,const ExpansionStorage<B> &b)
  {
    ExpansionStorage<bounded_capacity(A+B)> result;
    double sum[A+B];
    result.assign(::predicates_expansion_sum(a.size,a.data,b.size,b.data,sum),sum);
    return result;
  }
  template<int A,int B> inline ExpansionStorage<bounded_capacity(2*A*B)> operator*(
    const ExpansionStorage<A> &a,const ExpansionStorage<B> &b)
  {
    ExpansionStorage<bounded_capacity(2*A*B)> result;
    if(sign(a)==0 || sign(b)==0) return result;
    double scaled[2*A],sum[bounded_capacity(2*A*B)+2*A];
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
  template<typename Number> using Input=typename std::conditional<
    std::is_same<Number,Expansion>::value,ExpansionStorage<1>,Number>::type;
  template<typename Number> inline auto difference(const Point &a,const Point &b)
  {
    using Scalar=Input<Number>;
    using Difference=decltype(Scalar(a[0])-Scalar(b[0]));
    return Vector<Difference>{{Scalar(a[0])-Scalar(b[0]),Scalar(a[1])-Scalar(b[1]),Scalar(a[2])-Scalar(b[2])}};
  }
  template<> inline auto difference<Filter>(const Point &a,const Point &b)
  { return Vector<Filter>{{Filter(a[0]-b[0]),Filter(a[1]-b[1]),Filter(a[2]-b[2])}}; }
  template<typename A,typename B> inline auto dot(const Vector<A> &a,const Vector<B> &b)
  { return (a[0]*b[0]+a[1]*b[1])+a[2]*b[2]; }
  template<typename A,typename B> inline auto cross(const Vector<A> &a,const Vector<B> &b)
  {
    using ProductDifference=decltype(a[1]*b[2]-a[2]*b[1]);
    return Vector<ProductDifference>{{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}};
  }

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
  // The existing 256*epsilon bound exceeds gamma_30/(1-gamma_30), including
  // rounding of this multiplication. Reject underflowed/nonfinite error estimates.
  inline bool filter_error(const Filter &value,double &bound)
  {
    bound=256*std::numeric_limits<double>::epsilon()*value.permanent;
    return std::isfinite(value.value) &&
      (std::isnormal(bound) || (bound==0 && value.permanent==0));
  }
  inline bool distance_upper_bound(const Fraction<Filter> &fraction,double &upper)
  {
    double numerator_error,denominator_error;
    if(!filter_error(fraction.numerator,numerator_error) || !filter_error(fraction.denominator,denominator_error) ||
       fraction.denominator.value<=denominator_error ||
       (fraction.numerator.value<=numerator_error && fraction.numerator.permanent!=0)) return false;
    if(fraction.numerator.permanent==0) { upper=0; return true; }
    const double numerator=up(fraction.numerator.value+numerator_error);
    const double denominator=down(fraction.denominator.value-denominator_error);
    upper=up(numerator/denominator);
    return denominator>0 && std::isfinite(upper);
  }

  // Test the same vertex-first, edge-second regions as the public classifier,
  // reusing five dot products. -1 asks the existing classifier to refine.
  inline int filtered_feature_distance(const Point &q,const std::array<Point,3> &p,
    int &id,Fraction<Filter> &distance)
  {
    const auto r=difference<Filter>(q,p[0]);
    const auto e=difference<Filter>(p[1],p[0]);
    const auto f=difference<Filter>(p[2],p[0]);
    const auto re=dot(r,e),rf=dot(r,f);
    const auto vertex=[&](const Filter &first,const Filter &second)
    {
      const int a=sign(first),b=sign(second);
      if(a==1 || b==1) return 0;
      return a==2 || b==2 ? -1 : 1;
    };
    int selected=vertex(re,rf);
    if(selected<0) return -1;
    if(selected) { id=0; distance={dot(r,r),Filter(1)}; return 0; }
    const auto ee=dot(e,e),ef=dot(e,f);
    // Projections from b onto ba and bc, then from c onto ca and cb.
    const auto br=ee-re,bc=(rf-ef)+br;
    selected=vertex(br,bc);
    if(selected<0) return -1;
    if(selected)
    {
      const auto offset=difference<Filter>(q,p[1]);
      id=1; distance={dot(offset,offset),Filter(1)}; return 0;
    }
    const auto ff=dot(f,f),cr=ff-rf,cb=(re-ef)+cr;
    selected=vertex(cr,cb);
    if(selected<0) return -1;
    if(selected)
    {
      const auto offset=difference<Filter>(q,p[2]);
      id=2; distance={dot(offset,offset),Filter(1)}; return 0;
    }
    const auto edge=[&](const Filter &first,const Filter &last)
    {
      const int a=sign(first),b=sign(last);
      if(a==2 || b==2) return -1;
      return a>0 && b>0 ? 1 : 0;
    };
    selected=edge(re,br);
    if(selected<0) return -1;
    if(selected)
    {
      const int side=sign(ee*rf-re*ef);
      if(side==2) return -1;
      if(side<=0) { id=0; distance={dot(r,r)*ee-re*re,ee}; return 1; }
    }
    selected=edge(bc,cb);
    if(selected<0) return -1;
    if(selected)
    {
      const auto be=ee-ef,ce=ff-ef,length=be+ce;
      const int side=sign(length*br-bc*be);
      if(side==2) return -1;
      if(side<=0)
      {
        const auto offset=difference<Filter>(q,p[1]);
        id=1; distance={dot(offset,offset)*length-bc*bc,length}; return 1;
      }
    }
    selected=edge(cr,rf);
    if(selected<0) return -1;
    if(selected)
    {
      const int side=sign(ff*cb-cr*(ff-ef));
      if(side==2) return -1;
      if(side<=0)
      {
        const auto offset=difference<Filter>(q,p[2]);
        id=2; distance={dot(offset,offset)*ff-cr*cr,ff}; return 1;
      }
    }
    const auto normal=cross(e,f);
    const auto height=dot(r,normal);
    id=0; distance={height*height,dot(normal,normal)}; return 2;
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
  // Reuse the incumbent's rounded fraction and permanent. Refinement uses
  // fresh local interval/expansion storage and the original input coordinates.
  inline int compare(const Point &q,
    const std::array<Point,3> &a,int ad,int ai,const Fraction<Filter> &first,
    const std::array<Point,3> &b,int bd,int bi,const Fraction<Filter> &second)
  {
    if(ad==0 && bd==0 && a[ai]==b[bi]) return 0;
    if(ad==1 && bd==1 &&
       ((a[ai]==b[bi] && a[(ai+1)%3]==b[(bi+1)%3]) ||
        (a[ai]==b[(bi+1)%3] && a[(ai+1)%3]==b[bi]))) return 0;
    const int fast=sign(first.numerator*second.denominator-second.numerator*first.denominator);
    if(fast!=2) return fast;
    const int filtered=compare<Interval>(q,a,ad,ai,b,bd,bi);
    return filtered==2 ? compare<Expansion>(q,a,ad,ai,b,bd,bi) : filtered;
  }
  // Reuse exact feature classifications, computing fractions within this call.
  // All interval and expansion storage remains local to the predicate.
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

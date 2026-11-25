// Benchmark-only controls for certified closest-point filters.
#pragma once
#include "filtered_distance.h"
#include <igl/predicates/exactinit.h>
#include <vector>
namespace filter_profile
{
struct Counts
{
  std::size_t nodes=0,leaves=0,updates=0,feature_fallbacks=0,interval_bounds=0,
    filter_compares=0,enclosure_compares=0,interval_compares=0,expansion_compares=0;
};
template<bool FastUpper,bool CacheDistance,bool Fused,bool Enclosures,bool Audit=false,
  typename DerivedQ,typename DerivedV,typename DerivedF,typename DerivedB,typename DerivedL>
Eigen::Index search(
  const Eigen::MatrixBase<DerivedQ> &query,
  const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
  const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
  const Eigen::MatrixBase<DerivedL> &leaf,Counts &counts)
{
  static_assert(!Fused || CacheDistance);
  static_assert(!Enclosures || CacheDistance);
  igl::predicates::exactinit();
  const auto q=k::point(query);
  if(V.cols()!=3 || F.cols()!=3 || B1.cols()!=3 || B2.cols()!=3 ||
     B1.rows()!=leaf.size() || B2.rows()!=leaf.size())
    throw std::invalid_argument("Incompatible triangle tree dimensions");
  if(leaf.size()==0) return -1;
  const auto triangle=[&](Eigen::Index face)
  {
    if(face<0 || face>=F.rows()) throw std::out_of_range("Invalid triangle leaf id");
    std::array<k::Point,3> result;
    for(int i=0;i<3;++i)
    {
      const Eigen::Index vertex=F(face,i);
      if(vertex<0 || vertex>=V.rows()) throw std::out_of_range("Invalid triangle vertex index");
      result[i]=k::point(V.row(vertex));
    }
    return result;
  };
  const auto lower_bound=[&](Eigen::Index node)
  {
    double squared=0;
    for(int axis=0;axis<3;++axis)
    {
      const double lo=B1(node,axis),hi=B2(node,axis);
      double delta=0;
      if(q[axis]<lo) delta=lo-q[axis];
      else if(q[axis]>hi) delta=q[axis]-hi;
      squared+=delta*delta;
    }
    if(!std::isfinite(squared)) return 0.0;
    // At most five relative rounding factors per positive term. The slack
    // also covers rounding of this subtraction; four subnormal ulps cover
    // absolute errors when a squared term underflows. Never round a lower
    // bound upward into a pruning decision.
    return std::max(0.0,(squared-16*std::numeric_limits<double>::epsilon()*squared)
      -4*std::numeric_limits<double>::denorm_min());
  };
  const auto upper_bound=[&](const std::array<k::Point,3> &p,int dimension,int id)
  {
    double upper=std::numeric_limits<double>::infinity();
    for(const auto &corner:p)
    {
      const auto offset=k::difference<k::Interval>(q,corner);
      upper=std::min(upper,k::dot(offset,offset).hi);
    }
    const auto fraction=k::distance<k::Interval>(q,p,dimension,id);
    if(fraction.denominator.lo>0)
      upper=std::min(upper,k::up(fraction.numerator.hi/fraction.denominator.lo));
    return upper;
  };
  struct Node { Eigen::Index index; double lower; };
  std::vector<Node> active;
  active.reserve(64);
  active.push_back({0,0});
  Eigen::Index best=-1;
  std::array<k::Point,3> best_triangle;
  int best_dimension=0,best_id=0;
  Fraction best_distance;
  k::Interval best_bounds;
  bool best_bounded=false;
  double upper=std::numeric_limits<double>::infinity();
  while(!active.empty())
  {
    const Node current=active.back(); active.pop_back();
    if constexpr(Audit) ++counts.nodes;
    if(current.lower>upper || leaf(current.index)==-2) continue;
    if(leaf(current.index)>=0)
    {
      const Eigen::Index candidate=leaf(current.index);
      const auto candidate_triangle=triangle(candidate);
      if constexpr(Audit) ++counts.leaves;
      int candidate_id,candidate_dimension=-1;
      Fraction candidate_distance;
      k::Interval candidate_bounds;
      bool candidate_bounded=false;
      if constexpr(Fused)
        candidate_dimension=feature_distance(q,candidate_triangle,candidate_id,candidate_distance);
      if(candidate_dimension<0)
      {
        if constexpr(Audit && Fused) ++counts.feature_fallbacks;
        candidate_dimension=k::feature(q,candidate_triangle,candidate_id);
        if constexpr(CacheDistance)
          candidate_distance=k::distance<k::Filter>(q,candidate_triangle,candidate_dimension,candidate_id);
      }
      if constexpr(Enclosures)
        candidate_bounded=distance_bounds(candidate_distance,candidate_bounds);
      int order=-1;
      if(best>=0)
      {
        const auto &a=candidate_triangle,&b=best_triangle;
        const int ad=candidate_dimension,ai=candidate_id,bd=best_dimension,bi=best_id;
        if constexpr(!CacheDistance) order=k::compare(q,a,ad,ai,b,bd,bi);
        else
        {
          if((ad==0 && bd==0 && a[ai]==b[bi]) || (ad==1 && bd==1 &&
             ((a[ai]==b[bi] && a[(ai+1)%3]==b[(bi+1)%3]) ||
              (a[ai]==b[(bi+1)%3] && a[(ai+1)%3]==b[bi])))) order=0;
          else if(Enclosures && candidate_bounded && best_bounded &&
            (candidate_bounds.lo>best_bounds.hi || candidate_bounds.hi<best_bounds.lo))
          {
            order=candidate_bounds.lo>best_bounds.hi ? 1 : -1;
            if constexpr(Audit) ++counts.enclosure_compares;
          }
          else
          {
            if constexpr(Audit) ++counts.filter_compares;
            order=k::sign(candidate_distance.numerator*best_distance.denominator-
              best_distance.numerator*candidate_distance.denominator);
            if(order==2)
            {
              if constexpr(Audit) ++counts.interval_compares;
              order=k::compare<k::Interval>(q,a,ad,ai,b,bd,bi);
              if(order==2)
              {
                if constexpr(Audit) ++counts.expansion_compares;
                order=k::compare<k::Expansion>(q,a,ad,ai,b,bd,bi);
              }
            }
          }
        }
      }
      if(order<0 || (order==0 && candidate<best))
      {
        best=candidate; best_triangle=candidate_triangle;
        best_dimension=candidate_dimension; best_id=candidate_id;
        // An exact tie has the same distance, so the previous conservative
        // bound remains valid even when the smaller triangle id wins.
        if constexpr(CacheDistance) best_distance=candidate_distance;
        if constexpr(Enclosures) { best_bounds=candidate_bounds; best_bounded=candidate_bounded; }
        if(order<0)
        {
          bool bounded=false;
          if constexpr(Audit) ++counts.updates;
          if constexpr(Enclosures)
          {
            if(best_bounded) { upper=best_bounds.hi; bounded=true; }
          }
          else if constexpr(FastUpper)
          {
            if constexpr(!CacheDistance)
              candidate_distance=k::distance<k::Filter>(q,best_triangle,best_dimension,best_id);
            bounded=filter_profile::upper_bound(candidate_distance,upper);
          }
          if(!bounded)
          {
            upper=upper_bound(best_triangle,best_dimension,best_id);
            if constexpr(Audit) ++counts.interval_bounds;
          }
        }
      }
    }else
    {
      const Eigen::Index left=2*current.index+1,right=left+1;
      if(leaf(current.index)!=-1 || right>=leaf.size())
        throw std::invalid_argument("Invalid Eytzinger internal node");
      // Empty padding boxes need not contain initialized coordinates.
      const double left_lower=leaf(left)==-2 ? std::numeric_limits<double>::infinity() : lower_bound(left);
      const double right_lower=leaf(right)==-2 ? std::numeric_limits<double>::infinity() : lower_bound(right);
      if(left_lower<right_lower)
      { active.push_back({right,right_lower}); active.push_back({left,left_lower}); }
      else
      { active.push_back({left,left_lower}); active.push_back({right,right_lower}); }
    }
  }
  return best;
}

}

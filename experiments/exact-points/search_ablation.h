// Benchmark-only copy of exact_point_mesh_closest_triangle's traversal.
// The Rounded variants deliberately remove correctness guarantees. They are
// controls for cost attribution, not alternative public implementations.
#pragma once
#include <igl/predicates/exact_point_mesh_closest_triangle.h>
#include <igl/point_simplex_squared_distance.h>

namespace search_profile
{
namespace k=igl::predicates::exact_point_triangle_kernel;
struct Counts
{
  std::size_t nodes=0,leaves=0,updates=0,comparisons=0;
  std::size_t ties=0,shared_vertex_ties=0,shared_edge_ties=0,tie_updates=0;
  std::size_t compare_intervals=0,compare_expansions=0;
  std::size_t upper_intervals=0,upper_feature_intervals=0,upper_feature_expansions=0;
};
enum class FeatureCache { None, Coordinates, Topology };

template<typename Number> int classified_compare(const k::Point &q,
  const std::array<k::Point,3> &a,int ad,int ai,
  const std::array<k::Point,3> &b,int bd,int bi)
{
  const auto first=k::distance<Number>(q,a,ad,ai),second=k::distance<Number>(q,b,bd,bi);
  return k::sign(first.numerator*second.denominator-second.numerator*first.denominator);
}

// Keep the traversal, validation and conservative box lower bounds identical
// across all four policies. Collect counts in a separate, untimed instantiation.
template<bool ExactCompare,bool IntervalUpper,bool Audit=false,bool KeepTieUpper=false,
  FeatureCache Cache=FeatureCache::None>
Eigen::Index search(const Eigen::RowVector3d &query,
  const Eigen::MatrixXd &V,const Eigen::MatrixXi &F,
  const Eigen::MatrixXd &B1,const Eigen::MatrixXd &B2,
  const Eigen::VectorXi &leaf,Counts &counts)
{
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
    return std::max(0.0,(squared-16*std::numeric_limits<double>::epsilon()*squared)
      -4*std::numeric_limits<double>::denorm_min());
  };
  const auto upper_bound=[&](const std::array<k::Point,3> &p,int dimension,int id)
  {
    if constexpr(Audit) ++counts.upper_intervals;
    double upper=std::numeric_limits<double>::infinity();
    for(const auto &corner:p)
    {
      const auto offset=k::difference<k::Interval>(q,corner);
      upper=std::min(upper,k::dot(offset,offset).hi);
    }
    if constexpr(Cache==FeatureCache::None && Audit)
    {
      dimension=k::feature<k::Filter>(q,p,id);
      if(dimension<0)
      {
        ++counts.upper_feature_intervals;
        dimension=k::feature<k::Interval>(q,p,id);
        if(dimension<0)
        {
          ++counts.upper_feature_expansions;
          dimension=k::feature<k::Expansion>(q,p,id);
        }
      }
    }else if constexpr(Cache==FeatureCache::None) dimension=k::feature(q,p,id);
    const auto fraction=k::distance<k::Interval>(q,p,dimension,id);
    if(fraction.denominator.lo>0)
      upper=std::min(upper,k::up(fraction.numerator.hi/fraction.denominator.lo));
    return upper;
  };
  const auto compare=[&](const auto &candidate,const auto &best)
  {
    if constexpr(Audit)
    {
      ++counts.comparisons;
      const int fast=k::compare<k::Filter>(q,candidate,best);
      if(fast!=2) return fast;
      ++counts.compare_intervals;
      const int filtered=k::compare<k::Interval>(q,candidate,best);
      if(filtered!=2) return filtered;
      ++counts.compare_expansions;
      return k::compare<k::Expansion>(q,candidate,best);
    }else return k::compare(q,candidate,best);
  };
  struct Node { Eigen::Index index; double lower; };
  std::vector<Node> active;
  active.reserve(64);
  active.push_back({0,0});
  Eigen::Index best=-1;
  std::array<k::Point,3> best_triangle;
  int best_dimension=-1,best_id=-1;
  double upper=std::numeric_limits<double>::infinity();
  double best_squared=upper;
  while(!active.empty())
  {
    const Node current=active.back(); active.pop_back();
    if constexpr(Audit) ++counts.nodes;
    if(current.lower>upper || leaf(current.index)==-2) continue;
    if(leaf(current.index)>=0)
    {
      if constexpr(Audit) ++counts.leaves;
      const Eigen::Index candidate=leaf(current.index);
      const auto candidate_triangle=triangle(candidate);
      int candidate_dimension=-1,candidate_id=-1;
      if constexpr(Cache!=FeatureCache::None)
        candidate_dimension=k::feature(q,candidate_triangle,candidate_id);
      double squared=0;
      if constexpr(!ExactCompare)
      {
        Eigen::RowVector3d closest;
        igl::point_simplex_squared_distance<3>(query,V,F,candidate,squared,closest);
      }
      int order;
      if constexpr(Cache!=FeatureCache::None)
      {
        static_assert(ExactCompare && IntervalUpper,"Feature caching applies to the exact control");
        if(best<0) order=-1;
        else
        {
          if constexpr(Audit) ++counts.comparisons;
          const auto same_vertex=[&](int a,int b)
          {
            if constexpr(Cache==FeatureCache::Topology) return F(candidate,a)==F(best,b);
            else return candidate_triangle[a]==best_triangle[b];
          };
          const bool shared_vertex=candidate_dimension==0 && best_dimension==0 &&
            same_vertex(candidate_id,best_id);
          const bool shared_edge=candidate_dimension==1 && best_dimension==1 &&
            ((same_vertex(candidate_id,best_id) && same_vertex((candidate_id+1)%3,(best_id+1)%3)) ||
             (same_vertex(candidate_id,(best_id+1)%3) && same_vertex((candidate_id+1)%3,best_id)));
          if(shared_vertex || shared_edge) order=0;
          else
          {
            order=classified_compare<k::Filter>(q,candidate_triangle,candidate_dimension,candidate_id,
              best_triangle,best_dimension,best_id);
            if(order==2)
            {
              if constexpr(Audit) ++counts.compare_intervals;
              order=classified_compare<k::Interval>(q,candidate_triangle,candidate_dimension,candidate_id,
                best_triangle,best_dimension,best_id);
              if(order==2)
              {
                if constexpr(Audit) ++counts.compare_expansions;
                order=classified_compare<k::Expansion>(q,candidate_triangle,candidate_dimension,candidate_id,
                  best_triangle,best_dimension,best_id);
              }
            }
          }
        }
      }else if constexpr(ExactCompare)
        order=best<0 ? -1 : compare(candidate_triangle,best_triangle);
      else order=best<0 ? -1 : (squared>best_squared)-(squared<best_squared);
      if constexpr(Audit && ExactCompare)
        if(order==0)
        {
          ++counts.ties;
          // Classify again only in the untimed audit to count topological ties.
          int ai,bi;
          const int ad=k::feature(q,candidate_triangle,ai),bd=k::feature(q,best_triangle,bi);
          if(ad==0 && bd==0 && F(candidate,ai)==F(best,bi)) ++counts.shared_vertex_ties;
          if(ad==1 && bd==1 &&
             ((F(candidate,ai)==F(best,bi) && F(candidate,(ai+1)%3)==F(best,(bi+1)%3)) ||
              (F(candidate,ai)==F(best,(bi+1)%3) && F(candidate,(ai+1)%3)==F(best,bi))))
            ++counts.shared_edge_ties;
        }
      if(order<0 || (order==0 && candidate<best))
      {
        if constexpr(Audit) { ++counts.updates; counts.tie_updates+=order==0; }
        if constexpr(ExactCompare && !IntervalUpper)
        {
          Eigen::RowVector3d closest;
          igl::point_simplex_squared_distance<3>(query,V,F,candidate,squared,closest);
        }
        best=candidate; best_triangle=candidate_triangle; best_squared=squared;
        best_dimension=candidate_dimension; best_id=candidate_id;
        if constexpr(IntervalUpper)
        {
          if(!KeepTieUpper || order!=0) upper=upper_bound(best_triangle,best_dimension,best_id);
        }
        else upper=squared;
      }
    }else
    {
      const Eigen::Index left=2*current.index+1,right=left+1;
      if(leaf(current.index)!=-1 || right>=leaf.size())
        throw std::invalid_argument("Invalid Eytzinger internal node");
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

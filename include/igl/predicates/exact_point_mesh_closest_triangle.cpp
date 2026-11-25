// This file is part of libigl, a simple c++ geometry processing library.
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "exact_point_mesh_closest_triangle.h"
#include "exact_point_triangle_kernel.h"
#include "exactinit.h"
#include <vector>

template<typename DerivedQ,typename DerivedV,typename DerivedF,
  typename DerivedB,typename DerivedL>
IGL_INLINE Eigen::Index igl::predicates::exact_point_mesh_closest_triangle(
  const Eigen::MatrixBase<DerivedQ> &query,
  const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
  const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
  const Eigen::MatrixBase<DerivedL> &leaf)
{
  namespace k=exact_point_triangle_kernel;
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
  const auto upper_bound=[&](const std::array<k::Point,3> &p,int dimension,int id,
    const k::Fraction<k::Filter> &filtered)
  {
    double upper;
    if(k::distance_upper_bound(filtered,upper)) return upper;
    upper=std::numeric_limits<double>::infinity();
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
  k::Fraction<k::Filter> best_distance;
  double upper=std::numeric_limits<double>::infinity();
  while(!active.empty())
  {
    const Node current=active.back(); active.pop_back();
    if(current.lower>upper || leaf(current.index)==-2) continue;
    if(leaf(current.index)>=0)
    {
      const Eigen::Index candidate=leaf(current.index);
      const auto candidate_triangle=triangle(candidate);
      int candidate_id;
      k::Fraction<k::Filter> candidate_distance;
      int candidate_dimension=k::filtered_feature_distance(q,candidate_triangle,candidate_id,candidate_distance);
      if(candidate_dimension<0)
      {
        candidate_dimension=k::feature(q,candidate_triangle,candidate_id);
        candidate_distance=k::distance<k::Filter>(q,candidate_triangle,candidate_dimension,candidate_id);
      }
      const int order=best<0 ? -1 : k::compare(q,
        candidate_triangle,candidate_dimension,candidate_id,candidate_distance,
        best_triangle,best_dimension,best_id,best_distance);
      if(order<0 || (order==0 && candidate<best))
      {
        best=candidate; best_triangle=candidate_triangle;
        best_dimension=candidate_dimension; best_id=candidate_id;
        best_distance=candidate_distance;
        // An exact tie has the same distance, so the previous conservative
        // bound remains valid even when the smaller triangle id wins.
        if(order<0) upper=upper_bound(best_triangle,best_dimension,best_id,best_distance);
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

#ifdef IGL_STATIC_LIBRARY
template Eigen::Index igl::predicates::exact_point_mesh_closest_triangle<
  Eigen::RowVector3d,Eigen::MatrixXd,Eigen::MatrixXi,Eigen::MatrixXd,Eigen::VectorXi>(
  const Eigen::MatrixBase<Eigen::RowVector3d>&,const Eigen::MatrixBase<Eigen::MatrixXd>&,
  const Eigen::MatrixBase<Eigen::MatrixXi>&,const Eigen::MatrixBase<Eigen::MatrixXd>&,
  const Eigen::MatrixBase<Eigen::MatrixXd>&,const Eigen::MatrixBase<Eigen::VectorXi>&);
template Eigen::Index igl::predicates::exact_point_mesh_closest_triangle<
  Eigen::Vector3d,Eigen::MatrixXd,Eigen::MatrixXi,Eigen::MatrixXd,Eigen::VectorXi>(
  const Eigen::MatrixBase<Eigen::Vector3d>&,const Eigen::MatrixBase<Eigen::MatrixXd>&,
  const Eigen::MatrixBase<Eigen::MatrixXi>&,const Eigen::MatrixBase<Eigen::MatrixXd>&,
  const Eigen::MatrixBase<Eigen::MatrixXd>&,const Eigen::MatrixBase<Eigen::VectorXi>&);
#endif

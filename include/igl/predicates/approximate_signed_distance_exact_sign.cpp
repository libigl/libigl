// This file is part of libigl, a simple c++ geometry processing library.
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "approximate_signed_distance_exact_sign.h"
#include "point_mesh_closest_point_kernel.h"
#include "../parallel_for.h"
#include <exception>
#include <vector>

namespace igl { namespace predicates { namespace approximate_signed_distance_exact_sign_kernel
{
  template<typename V,typename F,typename B,typename L> inline void dimensions(
    const V &vertices,const F &faces,const B &lower,const B &upper,
    const L &leaves,const MeshFeatureExposure &exposure)
  {
    if(vertices.cols()!=3 || faces.cols()!=3 || lower.cols()!=3 || upper.cols()!=3 ||
       lower.rows()!=leaves.size() || upper.rows()!=leaves.size() ||
       exposure.edges.rows()!=faces.rows() || exposure.edges.cols()!=3 ||
       exposure.vertices.size()!=vertices.rows())
      throw std::invalid_argument("Incompatible signed distance tree/cache dimensions");
    if(faces.rows()==0 || leaves.size()==0)
      throw std::invalid_argument("An empty mesh/tree has no signed distance");
  }
  template<typename DerivedP,typename Evaluate> inline void batch(
    const Eigen::MatrixBase<DerivedP> &P,const Evaluate &evaluate,
    Eigen::VectorXd &S,Eigen::VectorXi &I,Eigen::MatrixXd &C,Eigen::MatrixXd &bary)
  {
    if(P.cols()!=3) throw std::invalid_argument("Signed distance queries must have three columns");
    S.resize(P.rows()); I.resize(P.rows()); C.resize(P.rows(),3); bary.resize(P.rows(),3);
    std::vector<std::exception_ptr> errors;
    igl::parallel_for(P.rows(),
      [&](const std::size_t threads) { errors.resize(threads); },
      [&](const Eigen::Index i,const std::size_t thread)
      {
        try
        {
          const Eigen::Matrix<typename DerivedP::Scalar,1,3> query=P.row(i);
          Eigen::Index triangle;
          Eigen::RowVector3d point,weights;
          S(i)=evaluate(query,triangle,point,weights);
          I(i)=static_cast<int>(triangle); C.row(i)=point; bary.row(i)=weights;
        }catch(...)
        {
          if(!errors[thread]) errors[thread]=std::current_exception();
        }
      },
      [&](const std::size_t thread)
      { if(errors[thread]) std::rethrow_exception(errors[thread]); },10000);
  }
}}}

template<typename DerivedQ,typename DerivedV,typename DerivedF,
  typename DerivedB,typename DerivedL>
IGL_INLINE double igl::predicates::approximate_signed_distance_exact_sign(
  const Eigen::MatrixBase<DerivedQ> &query,
  const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
  const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
  const Eigen::MatrixBase<DerivedL> &leaf,const MeshFeatureExposure &exposure,
  Eigen::Index &I,Eigen::RowVector3d &C,Eigen::RowVector3d &bary)
{
  namespace k=exact_point_triangle_kernel;
  namespace m=exact_point_mesh_kernel;
  approximate_signed_distance_exact_sign_kernel::dimensions(V,F,B1,B2,leaf,exposure);
  const auto closest=point_mesh_closest_point_kernel::closest_point(query,V,F,B1,B2,leaf);
  I=closest.triangle; C=closest.closest_point; bary=closest.barycentric_coordinates;
  const Eigen::Index face=closest.triangle;
  const int dimension=closest.dimension,id=closest.feature_id;
  const std::array<k::Point,3> p{{m::corner(V,F,face,0),
    m::corner(V,F,face,1),m::corner(V,F,face,2)}};
  const int label=dimension==0 ? exposure.vertices(F(face,id)) :
    dimension==1 ? exposure.edges(face,id) : 0;
  const int location=m::cached_location(k::point(query),p,dimension,id,label);
  return location==0 ? 0.0 :
    location*std::max(closest.distance,std::numeric_limits<double>::denorm_min());
}

template<typename DerivedQ,typename DerivedV,typename DerivedF,
  typename DerivedB,typename DerivedL,typename DerivedTT,typename DerivedTTi>
IGL_INLINE double igl::predicates::approximate_signed_distance_exact_sign(
  const Eigen::MatrixBase<DerivedQ> &query,
  const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
  const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
  const Eigen::MatrixBase<DerivedL> &leaf,const Eigen::MatrixBase<DerivedTT> &TT,
  const Eigen::MatrixBase<DerivedTTi> &TTi,
  Eigen::Index &I,Eigen::RowVector3d &C,Eigen::RowVector3d &bary)
{
  namespace k=exact_point_triangle_kernel;
  namespace m=exact_point_mesh_kernel;
  m::dimensions(V,F,TT,TTi);
  const auto closest=point_mesh_closest_point_kernel::closest_point(query,V,F,B1,B2,leaf);
  I=closest.triangle; C=closest.closest_point; bary=closest.barycentric_coordinates;
  const Eigen::Index face=closest.triangle;
  const std::array<k::Point,3> p{{m::corner(V,F,face,0),
    m::corner(V,F,face,1),m::corner(V,F,face,2)}};
  const int location=m::local_location(V,F,TT,TTi,k::point(query),p,
    face,closest.dimension,closest.feature_id);
  return location==0 ? 0.0 :
    location*std::max(closest.distance,std::numeric_limits<double>::denorm_min());
}

template<typename DerivedP,typename DerivedV,typename DerivedF,
  typename DerivedB,typename DerivedL>
IGL_INLINE void igl::predicates::approximate_signed_distance_exact_sign(
  const Eigen::MatrixBase<DerivedP> &P,
  const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
  const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
  const Eigen::MatrixBase<DerivedL> &leaf,const MeshFeatureExposure &exposure,
  Eigen::VectorXd &S,Eigen::VectorXi &I,Eigen::MatrixXd &C,Eigen::MatrixXd &bary)
{
  approximate_signed_distance_exact_sign_kernel::dimensions(V,F,B1,B2,leaf,exposure);
  approximate_signed_distance_exact_sign_kernel::batch(P,
    [&](const Eigen::Matrix<typename DerivedP::Scalar,1,3> &query,
      Eigen::Index &triangle,Eigen::RowVector3d &point,Eigen::RowVector3d &weights)
    { return approximate_signed_distance_exact_sign(query,V,F,B1,B2,leaf,exposure,triangle,point,weights); },
    S,I,C,bary);
}

template<typename DerivedP,typename DerivedV,typename DerivedF,
  typename DerivedB,typename DerivedL,typename DerivedTT,typename DerivedTTi>
IGL_INLINE void igl::predicates::approximate_signed_distance_exact_sign(
  const Eigen::MatrixBase<DerivedP> &P,
  const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
  const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
  const Eigen::MatrixBase<DerivedL> &leaf,const Eigen::MatrixBase<DerivedTT> &TT,
  const Eigen::MatrixBase<DerivedTTi> &TTi,Eigen::VectorXd &S,
  Eigen::VectorXi &I,Eigen::MatrixXd &C,Eigen::MatrixXd &bary)
{
  exact_point_mesh_kernel::dimensions(V,F,TT,TTi);
  if(B1.cols()!=3 || B2.cols()!=3 || B1.rows()!=leaf.size() || B2.rows()!=leaf.size())
    throw std::invalid_argument("Incompatible triangle tree dimensions");
  if(F.rows()==0 || leaf.size()==0)
    throw std::invalid_argument("An empty mesh/tree has no signed distance");
  approximate_signed_distance_exact_sign_kernel::batch(P,
    [&](const Eigen::Matrix<typename DerivedP::Scalar,1,3> &query,
      Eigen::Index &triangle,Eigen::RowVector3d &point,Eigen::RowVector3d &weights)
    { return approximate_signed_distance_exact_sign(query,V,F,B1,B2,leaf,TT,TTi,triangle,point,weights); },
    S,I,C,bary);
}

#ifdef IGL_STATIC_LIBRARY
#define IGL_SIGNED_DISTANCE_EXACT_SIGN(Q) \
template double igl::predicates::approximate_signed_distance_exact_sign< \
  Q,Eigen::MatrixXd,Eigen::MatrixXi,Eigen::MatrixXd,Eigen::VectorXi>( \
  const Eigen::MatrixBase<Q>&,const Eigen::MatrixBase<Eigen::MatrixXd>&, \
  const Eigen::MatrixBase<Eigen::MatrixXi>&,const Eigen::MatrixBase<Eigen::MatrixXd>&, \
  const Eigen::MatrixBase<Eigen::MatrixXd>&,const Eigen::MatrixBase<Eigen::VectorXi>&, \
  const igl::predicates::MeshFeatureExposure&,Eigen::Index&,Eigen::RowVector3d&,Eigen::RowVector3d&);
IGL_SIGNED_DISTANCE_EXACT_SIGN(Eigen::RowVector3d)
IGL_SIGNED_DISTANCE_EXACT_SIGN(Eigen::Vector3d)
#undef IGL_SIGNED_DISTANCE_EXACT_SIGN
template void igl::predicates::approximate_signed_distance_exact_sign<
  Eigen::MatrixXd,Eigen::MatrixXd,Eigen::MatrixXi,Eigen::MatrixXd,Eigen::VectorXi>(
  const Eigen::MatrixBase<Eigen::MatrixXd>&,const Eigen::MatrixBase<Eigen::MatrixXd>&,
  const Eigen::MatrixBase<Eigen::MatrixXi>&,const Eigen::MatrixBase<Eigen::MatrixXd>&,
  const Eigen::MatrixBase<Eigen::MatrixXd>&,const Eigen::MatrixBase<Eigen::VectorXi>&,
  const igl::predicates::MeshFeatureExposure&,Eigen::VectorXd&,
  Eigen::VectorXi&,Eigen::MatrixXd&,Eigen::MatrixXd&);
#define IGL_UNCACHED_SIGNED_DISTANCE(Q) \
template double igl::predicates::approximate_signed_distance_exact_sign< \
  Q,Eigen::MatrixXd,Eigen::MatrixXi,Eigen::MatrixXd,Eigen::VectorXi,Eigen::MatrixXi,Eigen::MatrixXi>( \
  const Eigen::MatrixBase<Q>&,const Eigen::MatrixBase<Eigen::MatrixXd>&, \
  const Eigen::MatrixBase<Eigen::MatrixXi>&,const Eigen::MatrixBase<Eigen::MatrixXd>&, \
  const Eigen::MatrixBase<Eigen::MatrixXd>&,const Eigen::MatrixBase<Eigen::VectorXi>&, \
  const Eigen::MatrixBase<Eigen::MatrixXi>&,const Eigen::MatrixBase<Eigen::MatrixXi>&, \
  Eigen::Index&,Eigen::RowVector3d&,Eigen::RowVector3d&);
IGL_UNCACHED_SIGNED_DISTANCE(Eigen::RowVector3d)
IGL_UNCACHED_SIGNED_DISTANCE(Eigen::Vector3d)
#undef IGL_UNCACHED_SIGNED_DISTANCE
template void igl::predicates::approximate_signed_distance_exact_sign<
  Eigen::MatrixXd,Eigen::MatrixXd,Eigen::MatrixXi,Eigen::MatrixXd,Eigen::VectorXi,Eigen::MatrixXi,Eigen::MatrixXi>(
  const Eigen::MatrixBase<Eigen::MatrixXd>&,const Eigen::MatrixBase<Eigen::MatrixXd>&,
  const Eigen::MatrixBase<Eigen::MatrixXi>&,const Eigen::MatrixBase<Eigen::MatrixXd>&,
  const Eigen::MatrixBase<Eigen::MatrixXd>&,const Eigen::MatrixBase<Eigen::VectorXi>&,
  const Eigen::MatrixBase<Eigen::MatrixXi>&,const Eigen::MatrixBase<Eigen::MatrixXi>&,
  Eigen::VectorXd&,Eigen::VectorXi&,Eigen::MatrixXd&,Eigen::MatrixXd&);
#endif

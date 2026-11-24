// This file is part of libigl, a simple c++ geometry processing library.
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "approximate_signed_distance_exact_sign.h"
#include "exact_point_mesh_kernel.h"
#include "../parallel_for.h"
#include <exception>

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
    Eigen::VectorXd &S,std::vector<PointMeshClosestPoint> &closest,Eigen::VectorXi &locations)
  {
    if(P.cols()!=3) throw std::invalid_argument("Signed distance queries must have three columns");
    S.resize(P.rows()); closest.resize(P.rows()); locations.resize(P.rows());
    std::vector<std::exception_ptr> errors;
    igl::parallel_for(P.rows(),
      [&](const std::size_t threads) { errors.resize(threads); },
      [&](const Eigen::Index i,const std::size_t thread)
      {
        try
        {
          const Eigen::Matrix<typename DerivedP::Scalar,1,3> query=P.row(i);
          PointMeshLocation location;
          S(i)=evaluate(query,closest[i],location);
          locations(i)=static_cast<int>(location);
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
  PointMeshClosestPoint &closest,PointMeshLocation &location)
{
  namespace k=exact_point_triangle_kernel;
  namespace m=exact_point_mesh_kernel;
  approximate_signed_distance_exact_sign_kernel::dimensions(V,F,B1,B2,leaf,exposure);
  closest=point_mesh_closest_point(query,V,F,B1,B2,leaf);
  const Eigen::Index face=closest.triangle;
  const int dimension=static_cast<int>(closest.feature),id=closest.feature_id;
  const std::array<k::Point,3> p{{m::corner(V,F,face,0),
    m::corner(V,F,face,1),m::corner(V,F,face,2)}};
  const int label=dimension==0 ? exposure.vertices(F(face,id)) :
    dimension==1 ? exposure.edges(face,id) : 0;
  location=static_cast<PointMeshLocation>(m::cached_location(k::point(query),p,dimension,id,label));
  return location==PointMeshLocation::ON_SURFACE ? 0.0 :
    static_cast<int>(location)*closest.distance;
}

template<typename DerivedQ,typename DerivedV,typename DerivedF,
  typename DerivedB,typename DerivedL,typename DerivedTT,typename DerivedTTi>
IGL_INLINE double igl::predicates::approximate_signed_distance_exact_sign(
  const Eigen::MatrixBase<DerivedQ> &query,
  const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
  const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
  const Eigen::MatrixBase<DerivedL> &leaf,const Eigen::MatrixBase<DerivedTT> &TT,
  const Eigen::MatrixBase<DerivedTTi> &TTi,
  PointMeshClosestPoint &closest,PointMeshLocation &location)
{
  namespace k=exact_point_triangle_kernel;
  namespace m=exact_point_mesh_kernel;
  m::dimensions(V,F,TT,TTi);
  closest=point_mesh_closest_point(query,V,F,B1,B2,leaf);
  const Eigen::Index face=closest.triangle;
  const std::array<k::Point,3> p{{m::corner(V,F,face,0),
    m::corner(V,F,face,1),m::corner(V,F,face,2)}};
  location=static_cast<PointMeshLocation>(m::local_location(V,F,TT,TTi,k::point(query),p,
    face,static_cast<int>(closest.feature),closest.feature_id));
  return location==PointMeshLocation::ON_SURFACE ? 0.0 :
    static_cast<int>(location)*closest.distance;
}

template<typename DerivedP,typename DerivedV,typename DerivedF,
  typename DerivedB,typename DerivedL>
IGL_INLINE void igl::predicates::approximate_signed_distance_exact_sign(
  const Eigen::MatrixBase<DerivedP> &P,
  const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
  const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
  const Eigen::MatrixBase<DerivedL> &leaf,const MeshFeatureExposure &exposure,
  Eigen::VectorXd &S,std::vector<PointMeshClosestPoint> &closest,Eigen::VectorXi &locations)
{
  approximate_signed_distance_exact_sign_kernel::dimensions(V,F,B1,B2,leaf,exposure);
  approximate_signed_distance_exact_sign_kernel::batch(P,
    [&](const Eigen::Matrix<typename DerivedP::Scalar,1,3> &query,
      PointMeshClosestPoint &point,PointMeshLocation &location)
    { return approximate_signed_distance_exact_sign(query,V,F,B1,B2,leaf,exposure,point,location); },
    S,closest,locations);
}

template<typename DerivedP,typename DerivedV,typename DerivedF,
  typename DerivedB,typename DerivedL,typename DerivedTT,typename DerivedTTi>
IGL_INLINE void igl::predicates::approximate_signed_distance_exact_sign(
  const Eigen::MatrixBase<DerivedP> &P,
  const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
  const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
  const Eigen::MatrixBase<DerivedL> &leaf,const Eigen::MatrixBase<DerivedTT> &TT,
  const Eigen::MatrixBase<DerivedTTi> &TTi,Eigen::VectorXd &S,
  std::vector<PointMeshClosestPoint> &closest,Eigen::VectorXi &locations)
{
  exact_point_mesh_kernel::dimensions(V,F,TT,TTi);
  if(B1.cols()!=3 || B2.cols()!=3 || B1.rows()!=leaf.size() || B2.rows()!=leaf.size())
    throw std::invalid_argument("Incompatible triangle tree dimensions");
  if(F.rows()==0 || leaf.size()==0)
    throw std::invalid_argument("An empty mesh/tree has no signed distance");
  approximate_signed_distance_exact_sign_kernel::batch(P,
    [&](const Eigen::Matrix<typename DerivedP::Scalar,1,3> &query,
      PointMeshClosestPoint &point,PointMeshLocation &location)
    { return approximate_signed_distance_exact_sign(query,V,F,B1,B2,leaf,TT,TTi,point,location); },
    S,closest,locations);
}

#ifdef IGL_STATIC_LIBRARY
#define IGL_SIGNED_DISTANCE_EXACT_SIGN(Q) \
template double igl::predicates::approximate_signed_distance_exact_sign< \
  Q,Eigen::MatrixXd,Eigen::MatrixXi,Eigen::MatrixXd,Eigen::VectorXi>( \
  const Eigen::MatrixBase<Q>&,const Eigen::MatrixBase<Eigen::MatrixXd>&, \
  const Eigen::MatrixBase<Eigen::MatrixXi>&,const Eigen::MatrixBase<Eigen::MatrixXd>&, \
  const Eigen::MatrixBase<Eigen::MatrixXd>&,const Eigen::MatrixBase<Eigen::VectorXi>&, \
  const igl::predicates::MeshFeatureExposure&,igl::predicates::PointMeshClosestPoint&, \
  igl::predicates::PointMeshLocation&);
IGL_SIGNED_DISTANCE_EXACT_SIGN(Eigen::RowVector3d)
IGL_SIGNED_DISTANCE_EXACT_SIGN(Eigen::Vector3d)
#undef IGL_SIGNED_DISTANCE_EXACT_SIGN
template void igl::predicates::approximate_signed_distance_exact_sign<
  Eigen::MatrixXd,Eigen::MatrixXd,Eigen::MatrixXi,Eigen::MatrixXd,Eigen::VectorXi>(
  const Eigen::MatrixBase<Eigen::MatrixXd>&,const Eigen::MatrixBase<Eigen::MatrixXd>&,
  const Eigen::MatrixBase<Eigen::MatrixXi>&,const Eigen::MatrixBase<Eigen::MatrixXd>&,
  const Eigen::MatrixBase<Eigen::MatrixXd>&,const Eigen::MatrixBase<Eigen::VectorXi>&,
  const igl::predicates::MeshFeatureExposure&,Eigen::VectorXd&,
  std::vector<igl::predicates::PointMeshClosestPoint>&,Eigen::VectorXi&);
#define IGL_UNCACHED_SIGNED_DISTANCE(Q) \
template double igl::predicates::approximate_signed_distance_exact_sign< \
  Q,Eigen::MatrixXd,Eigen::MatrixXi,Eigen::MatrixXd,Eigen::VectorXi,Eigen::MatrixXi,Eigen::MatrixXi>( \
  const Eigen::MatrixBase<Q>&,const Eigen::MatrixBase<Eigen::MatrixXd>&, \
  const Eigen::MatrixBase<Eigen::MatrixXi>&,const Eigen::MatrixBase<Eigen::MatrixXd>&, \
  const Eigen::MatrixBase<Eigen::MatrixXd>&,const Eigen::MatrixBase<Eigen::VectorXi>&, \
  const Eigen::MatrixBase<Eigen::MatrixXi>&,const Eigen::MatrixBase<Eigen::MatrixXi>&, \
  igl::predicates::PointMeshClosestPoint&,igl::predicates::PointMeshLocation&);
IGL_UNCACHED_SIGNED_DISTANCE(Eigen::RowVector3d)
IGL_UNCACHED_SIGNED_DISTANCE(Eigen::Vector3d)
#undef IGL_UNCACHED_SIGNED_DISTANCE
template void igl::predicates::approximate_signed_distance_exact_sign<
  Eigen::MatrixXd,Eigen::MatrixXd,Eigen::MatrixXi,Eigen::MatrixXd,Eigen::VectorXi,Eigen::MatrixXi,Eigen::MatrixXi>(
  const Eigen::MatrixBase<Eigen::MatrixXd>&,const Eigen::MatrixBase<Eigen::MatrixXd>&,
  const Eigen::MatrixBase<Eigen::MatrixXi>&,const Eigen::MatrixBase<Eigen::MatrixXd>&,
  const Eigen::MatrixBase<Eigen::MatrixXd>&,const Eigen::MatrixBase<Eigen::VectorXi>&,
  const Eigen::MatrixBase<Eigen::MatrixXi>&,const Eigen::MatrixBase<Eigen::MatrixXi>&,
  Eigen::VectorXd&,std::vector<igl::predicates::PointMeshClosestPoint>&,Eigen::VectorXi&);
#endif

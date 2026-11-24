// This file is part of libigl, a simple c++ geometry processing library.
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "approximate_signed_distance_exact_sign.h"
#include "exact_point_mesh_kernel.h"

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

template<typename DerivedP,typename DerivedV,typename DerivedF,
  typename DerivedB,typename DerivedL>
IGL_INLINE void igl::predicates::approximate_signed_distance_exact_sign(
  const Eigen::MatrixBase<DerivedP> &P,
  const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
  const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
  const Eigen::MatrixBase<DerivedL> &leaf,const MeshFeatureExposure &exposure,
  Eigen::VectorXd &S,std::vector<PointMeshClosestPoint> &closest,Eigen::VectorXi &locations)
{
  if(P.cols()!=3) throw std::invalid_argument("Signed distance queries must have three columns");
  approximate_signed_distance_exact_sign_kernel::dimensions(V,F,B1,B2,leaf,exposure);
  S.resize(P.rows()); closest.resize(P.rows()); locations.resize(P.rows());
  for(Eigen::Index i=0;i<P.rows();++i)
  {
    // Preserve the input scalar type; the single-query predicate validates it.
    const Eigen::Matrix<typename DerivedP::Scalar,1,3> query=P.row(i);
    PointMeshLocation location;
    S(i)=approximate_signed_distance_exact_sign(query,V,F,B1,B2,leaf,exposure,closest[i],location);
    locations(i)=static_cast<int>(location);
  }
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
#endif

// This file is part of libigl, a simple c++ geometry processing library.
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "exact_point_mesh_classify.h"
#include "exact_point_mesh_kernel.h"
#include "exactinit.h"

template<typename DerivedQ,typename DerivedV,typename DerivedF,
  typename DerivedTT,typename DerivedTTi>
IGL_INLINE igl::predicates::PointMeshLocation igl::predicates::exact_point_mesh_classify(
  const Eigen::MatrixBase<DerivedQ> &query,
  const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
  Eigen::Index closest_triangle,const Eigen::MatrixBase<DerivedTT> &TT,
  const Eigen::MatrixBase<DerivedTTi> &TTi)
{
  namespace m=exact_point_mesh_kernel;
  namespace k=exact_point_triangle_kernel;
  igl::predicates::exactinit();
  m::dimensions(V,F,TT,TTi);
  const auto q=k::point(query);
  const std::array<k::Point,3> p{{m::corner(V,F,closest_triangle,0),
    m::corner(V,F,closest_triangle,1),m::corner(V,F,closest_triangle,2)}};
  int id;
  const int dimension=k::feature(q,p,id);
  if(m::on_surface(q,p,dimension,id)) return PointMeshLocation::ON_SURFACE;
  if(dimension==2) return static_cast<PointMeshLocation>(-k::orientation(p[0],p[1],p[2],q));
  m::Corner edge{closest_triangle,id};
  if(dimension==0)
  {
    const auto ring=m::one_ring(F,TT,TTi,edge);
    for(const auto &candidate:ring)
    {
      const auto a=m::corner(V,F,candidate.face,(candidate.corner+1)%3);
      const auto b=m::corner(V,F,edge.face,(edge.corner+1)%3);
      if(k::perpendicular_compare(q,p[id],a,b)<0) edge=candidate;
    }
  }
  const int fold=m::fold(V,F,TT,TTi,edge);
  if(fold!=0) return static_cast<PointMeshLocation>(fold);
  return static_cast<PointMeshLocation>(-k::orientation(m::corner(V,F,edge.face,0),
    m::corner(V,F,edge.face,1),m::corner(V,F,edge.face,2),q));
}

template<typename DerivedQ,typename DerivedV,typename DerivedF>
IGL_INLINE igl::predicates::PointMeshLocation igl::predicates::exact_point_mesh_classify(
  const Eigen::MatrixBase<DerivedQ> &query,
  const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
  Eigen::Index closest_triangle,const MeshFeatureExposure &exposure)
{
  namespace m=exact_point_mesh_kernel;
  namespace k=exact_point_triangle_kernel;
  igl::predicates::exactinit();
  if(V.cols()!=3 || F.cols()!=3 || exposure.edges.rows()!=F.rows() ||
     exposure.edges.cols()!=3 || exposure.vertices.size()!=V.rows())
    throw std::invalid_argument("Incompatible feature exposure cache dimensions");
  const auto q=k::point(query);
  const std::array<k::Point,3> p{{m::corner(V,F,closest_triangle,0),
    m::corner(V,F,closest_triangle,1),m::corner(V,F,closest_triangle,2)}};
  int id;
  const int dimension=k::feature(q,p,id);
  const int exposure_label=dimension==0 ? exposure.vertices(F(closest_triangle,id)) :
    dimension==1 ? exposure.edges(closest_triangle,id) : 0;
  return static_cast<PointMeshLocation>(m::cached_location(q,p,dimension,id,exposure_label));
}

#ifdef IGL_STATIC_LIBRARY
#define IGL_POINT_MESH_CLASSIFY(Q) \
template igl::predicates::PointMeshLocation igl::predicates::exact_point_mesh_classify< \
  Q,Eigen::MatrixXd,Eigen::MatrixXi,Eigen::MatrixXi,Eigen::MatrixXi>( \
  const Eigen::MatrixBase<Q>&,const Eigen::MatrixBase<Eigen::MatrixXd>&, \
  const Eigen::MatrixBase<Eigen::MatrixXi>&,Eigen::Index,const Eigen::MatrixBase<Eigen::MatrixXi>&, \
  const Eigen::MatrixBase<Eigen::MatrixXi>&); \
template igl::predicates::PointMeshLocation igl::predicates::exact_point_mesh_classify< \
  Q,Eigen::MatrixXd,Eigen::MatrixXi>(const Eigen::MatrixBase<Q>&, \
  const Eigen::MatrixBase<Eigen::MatrixXd>&,const Eigen::MatrixBase<Eigen::MatrixXi>&, \
  Eigen::Index,const igl::predicates::MeshFeatureExposure&);
IGL_POINT_MESH_CLASSIFY(Eigen::RowVector3d)
IGL_POINT_MESH_CLASSIFY(Eigen::Vector3d)
#undef IGL_POINT_MESH_CLASSIFY
#endif

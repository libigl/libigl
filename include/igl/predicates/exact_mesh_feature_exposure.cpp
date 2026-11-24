// This file is part of libigl, a simple c++ geometry processing library.
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "exact_mesh_feature_exposure.h"
#include "exact_point_mesh_kernel.h"
#include "exactinit.h"

template<typename DerivedV,typename DerivedF,typename DerivedTT,typename DerivedTTi>
IGL_INLINE void igl::predicates::exact_mesh_feature_exposure(
  const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
  const Eigen::MatrixBase<DerivedTT> &TT,const Eigen::MatrixBase<DerivedTTi> &TTi,
  MeshFeatureExposure &exposure)
{
  namespace m=exact_point_mesh_kernel;
  namespace k=exact_point_triangle_kernel;
  igl::predicates::exactinit();
  m::dimensions(V,F,TT,TTi);
  exposure.edges.resize(F.rows(),3);
  exposure.vertices.setConstant(V.rows(),int(FeatureExposure::UNREACHABLE));
  std::vector<m::Corner> seeds(V.rows(),{-1,0});
  std::vector<Eigen::Index> incidence(V.rows(),0);
  for(Eigen::Index face=0;face<F.rows();++face)
  {
    const auto a=m::corner(V,F,face,0),b=m::corner(V,F,face,1),c=m::corner(V,F,face,2);
    // Reject degenerate faces before using a fold or a tangent plane.
    const auto normal=k::cross(k::difference<k::Expansion>(b,a),k::difference<k::Expansion>(c,a));
    if(k::sign(normal[0])==0 && k::sign(normal[1])==0 && k::sign(normal[2])==0)
      throw std::invalid_argument("Feature exposure requires nondegenerate faces");
    for(int i=0;i<3;++i)
    {
      seeds[F(face,i)]={face,i}; ++incidence[F(face,i)];
      exposure.edges(face,i)=m::fold(V,F,TT,TTi,{face,i});
    }
  }
  for(Eigen::Index vertex=0;vertex<V.rows();++vertex)
  {
    if(seeds[vertex].face<0) continue;
    const auto ring=m::one_ring(F,TT,TTi,seeds[vertex]);
    if(static_cast<Eigen::Index>(ring.size())!=incidence[vertex])
      throw std::invalid_argument("Feature exposure requires a single one-ring per vertex");
    const auto v=k::point(V.row(vertex));
    std::vector<k::Point> neighbors;
    for(const auto &entry:ring) neighbors.push_back(m::corner(V,F,entry.face,(entry.corner+1)%3));
    const auto seed=ring.front();
    const auto a=m::corner(V,F,seed.face,(seed.corner+1)%3);
    const auto b=m::corner(V,F,seed.face,(seed.corner+2)%3);
    bool planar=true;
    for(const auto &p:neighbors) if(k::orientation(v,a,b,p)!=0) { planar=false; break; }
    if(planar) { exposure.vertices(vertex)=int(FeatureExposure::FLAT); continue; }
    bool found=false;
    for(std::size_t i=0;i<ring.size() && !found;++i)
    {
      for(std::size_t j=i+1;j<ring.size() && !found;++j)
      {
        int side=0;
        bool supporting=true;
        for(const auto &p:neighbors)
        {
          const int orientation=k::orientation(v,neighbors[i],neighbors[j],p);
          if(orientation!=0)
          {
            if(side!=0 && side!=orientation) { supporting=false; break; }
            side=orientation;
          }
        }
        // An all-zero test is a collinear pair in this nonplanar one-ring.
        if(!supporting || side==0) continue;
        const auto edge=ring[i];
        int label=exposure.edges(edge.face,edge.corner);
        if(label==0)
        {
          const auto third=m::corner(V,F,edge.face,(edge.corner+2)%3);
          for(const auto &p:neighbors)
          {
            label=k::orientation(v,neighbors[i],third,p);
            if(label!=0) break;
          }
          if(label==0) throw std::logic_error("Nonplanar one-ring has no point off a supporting face");
        }
        exposure.vertices(vertex)=label;
        found=true;
      }
    }
  }
}
#ifdef IGL_STATIC_LIBRARY
template void igl::predicates::exact_mesh_feature_exposure<
  Eigen::MatrixXd,Eigen::MatrixXi,Eigen::MatrixXi,Eigen::MatrixXi>(
  const Eigen::MatrixBase<Eigen::MatrixXd>&,const Eigen::MatrixBase<Eigen::MatrixXi>&,
  const Eigen::MatrixBase<Eigen::MatrixXi>&,const Eigen::MatrixBase<Eigen::MatrixXi>&,
  igl::predicates::MeshFeatureExposure&);
#endif

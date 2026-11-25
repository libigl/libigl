// This file is part of libigl, a simple c++ geometry processing library.
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PREDICATES_EXACT_POINT_MESH_KERNEL_H
#define IGL_PREDICATES_EXACT_POINT_MESH_KERNEL_H
#include "exact_point_triangle_kernel.h"
#include <vector>
namespace igl { namespace predicates { namespace exact_point_mesh_kernel
{
  namespace k=exact_point_triangle_kernel;
  struct Corner { Eigen::Index face; int corner; };
  template<typename V,typename F> inline k::Point corner(const V &vertices,const F &faces,
    Eigen::Index face,int index)
  {
    if(face<0 || face>=faces.rows()) throw std::out_of_range("Invalid closest triangle id");
    const Eigen::Index vertex=faces(face,index);
    if(vertex<0 || vertex>=vertices.rows()) throw std::out_of_range("Invalid triangle vertex index");
    return k::point(vertices.row(vertex));
  }
  template<typename F,typename TT,typename TTi> inline Corner adjacent(
    const F &faces,const TT &neighbors,const TTi &neighbor_edges,Corner edge)
  {
    const Eigen::Index face=neighbors(edge.face,edge.corner);
    const int index=neighbor_edges(edge.face,edge.corner);
    if(face<0 || face>=faces.rows() || index<0 || index>=3 ||
       neighbors(face,index)!=edge.face || neighbor_edges(face,index)!=edge.corner ||
       faces(edge.face,edge.corner)!=faces(face,(index+1)%3) ||
       faces(edge.face,(edge.corner+1)%3)!=faces(face,index))
      throw std::invalid_argument("Expected closed, consistently oriented manifold adjacency");
    return {face,index};
  }
  template<typename F,typename TT,typename TTi> inline std::vector<Corner> one_ring(
    const F &faces,const TT &neighbors,const TTi &neighbor_edges,Corner seed)
  {
    std::vector<Corner> ring;
    Corner current=seed;
    do
    {
      if(ring.size()>=static_cast<std::size_t>(faces.rows()))
        throw std::invalid_argument("Vertex adjacency does not form a single cycle");
      ring.push_back(current);
      const Corner other=adjacent(faces,neighbors,neighbor_edges,current);
      current={other.face,(other.corner+1)%3};
    }while(current.face!=seed.face || current.corner!=seed.corner);
    return ring;
  }
  template<typename V,typename F,typename TT,typename TTi> inline void dimensions(
    const V &vertices,const F &faces,const TT &neighbors,const TTi &neighbor_edges)
  {
    if(vertices.cols()!=3 || faces.cols()!=3 || neighbors.cols()!=3 ||
       neighbor_edges.cols()!=3 || neighbors.rows()!=faces.rows() || neighbor_edges.rows()!=faces.rows())
      throw std::invalid_argument("Incompatible mesh adjacency dimensions");
  }
  template<typename V,typename F,typename TT,typename TTi> inline int fold(
    const V &vertices,const F &faces,const TT &neighbors,const TTi &neighbor_edges,Corner edge)
  {
    const Corner other=adjacent(faces,neighbors,neighbor_edges,edge);
    const auto a=corner(vertices,faces,edge.face,edge.corner),
      b=corner(vertices,faces,edge.face,(edge.corner+1)%3),
      c=corner(vertices,faces,edge.face,(edge.corner+2)%3),
      d=corner(vertices,faces,other.face,(other.corner+2)%3);
    const int label=k::orientation(a,b,c,d);
    // A zero fold can also mean a degenerate incident face. Only this branch
    // needs area checks before interpreting the fold as a flat surface.
    if(label==0 && (k::degenerate(a,b,c) || k::degenerate(a,b,d)))
      throw std::invalid_argument("Mesh signing requires nondegenerate faces");
    return label;
  }
  inline bool on_surface(const k::Point &q,const std::array<k::Point,3> &p,int dimension,int id)
  {
    if(dimension==0) return q==p[id];
    if(dimension==2) return k::orientation(p[0],p[1],p[2],q)==0;
    const int fast=k::sign(k::distance<k::Filter>(q,p,dimension,id).numerator);
    if(fast!=2) return fast==0;
    const int filtered=k::sign(k::distance<k::Interval>(q,p,dimension,id).numerator);
    if(filtered!=2) return filtered==0;
    std::array<k::Point,4> scaled{{q,p[0],p[1],p[2]}};
    k::rescale<6>(scaled);
    return k::sign(k::distance<k::Expansion>(scaled[0],{{scaled[1],scaled[2],scaled[3]}},dimension,id).numerator)==0;
  }
  template<typename V,typename F,typename TT,typename TTi> inline int local_location(
    const V &vertices,const F &faces,const TT &neighbors,const TTi &neighbor_edges,
    const k::Point &q,const std::array<k::Point,3> &p,Eigen::Index closest_triangle,
    int dimension,int id)
  {
    if(on_surface(q,p,dimension,id)) return 0;
    if(dimension==2) return -k::orientation(p[0],p[1],p[2],q);
    Corner edge{closest_triangle,id};
    if(dimension==0)
    {
      const auto ring=one_ring(faces,neighbors,neighbor_edges,edge);
      // The first edge is already the incumbent. Its self-comparison is an
      // exact zero and would force unnecessary expansion arithmetic.
      for(std::size_t i=1;i<ring.size();++i)
      {
        const auto &candidate=ring[i];
        const auto a=corner(vertices,faces,candidate.face,(candidate.corner+1)%3);
        const auto b=corner(vertices,faces,edge.face,(edge.corner+1)%3);
        if(k::perpendicular_compare(q,p[id],a,b)<0) edge=candidate;
      }
    }
    const int fold_label=fold(vertices,faces,neighbors,neighbor_edges,edge);
    if(fold_label!=0) return fold_label;
    const int side=-k::orientation(corner(vertices,faces,edge.face,0),
      corner(vertices,faces,edge.face,1),corner(vertices,faces,edge.face,2),q);
    if(side==0) throw std::invalid_argument("Noncoincident closest feature has no oriented side");
    return side;
  }
  inline int cached_location(const k::Point &q,const std::array<k::Point,3> &p,
    int dimension,int id,int exposure_label)
  {
    if(on_surface(q,p,dimension,id)) return 0;
    if(exposure_label==2)
      throw std::invalid_argument("A noncoincident query cannot have an unreachable closest vertex");
    if(exposure_label!=0 && exposure_label!=-1 && exposure_label!=1)
      throw std::invalid_argument("Invalid feature exposure label");
    if(exposure_label!=0) return exposure_label;
    return -k::orientation(p[0],p[1],p[2],q);
  }
}}}
#endif

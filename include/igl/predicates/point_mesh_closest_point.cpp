// This file is part of libigl, a simple c++ geometry processing library.
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "point_mesh_closest_point.h"
#include "exact_point_mesh_closest_triangle.h"
#include "exact_point_mesh_kernel.h"

namespace igl { namespace predicates { namespace point_mesh_closest_point_kernel
{
  namespace k=exact_point_triangle_kernel;
  template<typename Number> struct Construction
  {
    Number numerator,denominator;
    std::array<Number,3> barycentric,point;
  };
  template<typename Number> inline Construction<Number> construct(
    const k::Point &q,const std::array<k::Point,3> &p,int dimension,int id)
  {
    Construction<Number> result;
    const auto r=k::difference<Number>(q,p[id]);
    if(dimension==0)
    {
      result.numerator=k::dot(r,r); result.denominator=Number(1);
      result.barycentric[id]=Number(1);
    }else
    {
      const auto e=k::difference<Number>(p[(id+1)%3],p[id]);
      if(dimension==1)
      {
        result.denominator=k::dot(e,e);
        result.barycentric[(id+1)%3]=k::dot(r,e);
        // Cross products avoid subtracting two large squared distances for
        // a query very close to an edge. These are the same exact fraction.
        const auto normal=k::cross(r,e);
        result.numerator=k::dot(normal,normal);
      }else
      {
        const auto f=k::difference<Number>(p[(id+2)%3],p[id]);
        const auto normal=k::cross(e,f);
        result.denominator=k::dot(normal,normal);
        const auto height=k::dot(r,normal);
        result.numerator=height*height;
        result.barycentric[(id+1)%3]=k::dot(k::cross(r,f),normal);
        result.barycentric[(id+2)%3]=k::dot(k::cross(e,r),normal);
      }
      result.barycentric[id]=result.denominator-
        result.barycentric[(id+1)%3]-result.barycentric[(id+2)%3];
    }
    for(int axis=0;axis<3;++axis)
      result.point[axis]=(k::Input<Number>(p[0][axis])*result.barycentric[0]+
        k::Input<Number>(p[1][axis])*result.barycentric[1])+k::Input<Number>(p[2][axis])*result.barycentric[2];
    return result;
  }
  inline bool well_conditioned(const k::Filter &value)
  {
    // The permanent bound used by the metric predicates also covers these
    // degree-six constructions. A 1024:1 cancellation limit bounds the
    // accepted component's relative error by about 6e-11 (256*epsilon*1024),
    // before division. This is a numerical approximation policy, separate
    // from all exact feature and sign decisions. Smaller thresholds trigger
    // costly expansions for ordinary projections far from small triangles.
    return std::isfinite(value.value) && std::isfinite(value.permanent) &&
      std::abs(value.value)>=value.permanent/1024;
  }
  inline double estimate(const k::Filter &value) { return value.value; }
  inline double estimate(const k::Expansion &value)
  {
    double sum=0;
    for(int i=0;i<value.size;++i) sum+=value.data[i];
    return sum;
  }
  template<typename Number> inline void approximate(
    const Construction<Number> &construction,PointMeshClosestPoint &result)
  {
    const double numerator=estimate(construction.numerator);
    const double denominator=estimate(construction.denominator);
    if(!(denominator>0) || !std::isfinite(denominator) ||
       !(numerator>=0) || !std::isfinite(numerator))
      throw std::overflow_error("Closest point construction exceeded its arithmetic range");
    result.squared_distance=numerator/denominator;
    // Do not square or divide first: the squared result may underflow even
    // when both the fraction and the distance are representable.
    result.distance=std::sqrt(numerator)/std::sqrt(denominator);
    for(int i=0;i<3;++i)
    {
      result.barycentric_coordinates(i)=estimate(construction.barycentric[i])/denominator;
      result.closest_point(i)=estimate(construction.point[i])/denominator;
    }
  }
  // Keep the large return slot off the filtered query's frame, including in
  // header-only builds. Caller-owned 512 KiB worker stacks must also suffice
  // when a construction needs exact arithmetic.
#if defined(_MSC_VER)
  __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
  __attribute__((noinline))
#endif
  inline void construct_exact(const k::Point &q,const std::array<k::Point,3> &p,
    int dimension,int id,PointMeshClosestPoint &result)
  {
    std::array<k::Point,4> scaled{{q,p[0],p[1],p[2]}};
    const int shift=k::rescale<6>(scaled);
    approximate(construct<k::Expansion>(scaled[0],{{scaled[1],scaled[2],scaled[3]}},dimension,id),result);
    result.distance=std::ldexp(result.distance,-shift);
    result.squared_distance=std::ldexp(result.squared_distance,-2*shift);
    for(int axis=0;axis<3;++axis) result.closest_point(axis)=std::ldexp(result.closest_point(axis),-shift);
  }
}}}

template<typename DerivedQ,typename DerivedV,typename DerivedF,
  typename DerivedB,typename DerivedL>
IGL_INLINE igl::predicates::PointMeshClosestPoint igl::predicates::point_mesh_closest_point(
  const Eigen::MatrixBase<DerivedQ> &query,
  const Eigen::MatrixBase<DerivedV> &V,const Eigen::MatrixBase<DerivedF> &F,
  const Eigen::MatrixBase<DerivedB> &B1,const Eigen::MatrixBase<DerivedB> &B2,
  const Eigen::MatrixBase<DerivedL> &leaf)
{
  namespace k=exact_point_triangle_kernel;
  namespace m=exact_point_mesh_kernel;
  namespace c=point_mesh_closest_point_kernel;
  PointMeshClosestPoint result;
  result.triangle=exact_point_mesh_closest_triangle(query,V,F,B1,B2,leaf);
  if(result.triangle<0) throw std::invalid_argument("An empty mesh/tree has no closest point");
  const auto q=k::point(query);
  const std::array<k::Point,3> p{{m::corner(V,F,result.triangle,0),
    m::corner(V,F,result.triangle,1),m::corner(V,F,result.triangle,2)}};
  const int dimension=k::feature(q,p,result.feature_id);
  result.feature=static_cast<ClosestPointFeature>(dimension);
  const auto filtered=c::construct<k::Filter>(q,p,dimension,result.feature_id);
  bool reliable=c::well_conditioned(filtered.numerator) && c::well_conditioned(filtered.denominator);
  for(int i=0;i<3;++i)
    reliable=reliable && c::well_conditioned(filtered.barycentric[i]) && c::well_conditioned(filtered.point[i]);
  if(reliable) c::approximate(filtered,result);
  else c::construct_exact(q,p,dimension,result.feature_id,result);
  return result;
}

#ifdef IGL_STATIC_LIBRARY
#define IGL_POINT_MESH_CLOSEST_POINT(Q) \
template igl::predicates::PointMeshClosestPoint igl::predicates::point_mesh_closest_point< \
  Q,Eigen::MatrixXd,Eigen::MatrixXi,Eigen::MatrixXd,Eigen::VectorXi>( \
  const Eigen::MatrixBase<Q>&,const Eigen::MatrixBase<Eigen::MatrixXd>&, \
  const Eigen::MatrixBase<Eigen::MatrixXi>&,const Eigen::MatrixBase<Eigen::MatrixXd>&, \
  const Eigen::MatrixBase<Eigen::MatrixXd>&,const Eigen::MatrixBase<Eigen::VectorXi>&);
IGL_POINT_MESH_CLOSEST_POINT(Eigen::RowVector3d)
IGL_POINT_MESH_CLOSEST_POINT(Eigen::Vector3d)
#undef IGL_POINT_MESH_CLOSEST_POINT
#endif

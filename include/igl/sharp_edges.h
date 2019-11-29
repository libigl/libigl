#ifndef IGL_SHARP_EDGES_H
#define IGL_SHARP_EDGES_H

#include <igl/igl_inline.h>
#include <Eigen/Core>

namespace igl
{
  // SHARP_EDGES Given a mesh, compute sharp edges.
  //
  // Inputs:
  //   V  #V by 3 list of vertex positions
  //   F  #F by 3 list of triangle mesh indices into V
  //   angle  dihedral angle considered to sharp (e.g., igl::PI * 0.11)
  // Outputs:
  //   SE  #SE by 2 list of edge indices into V
  //
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedSE,
    typename DerivedE,
    typename DeriveduE,
    typename DerivedEMAP,
    typename uE2Etype,
    typename sharptype>
  IGL_INLINE void sharp_edges(
    const Eigen::MatrixBase<DerivedV> & V,
    const Eigen::MatrixBase<DerivedF> & F,
    const typename DerivedV::Scalar angle,
    Eigen::PlainObjectBase<DerivedSE> & SE,
    Eigen::PlainObjectBase<DerivedE> & E,
    Eigen::PlainObjectBase<DeriveduE> & uE,
    Eigen::PlainObjectBase<DerivedEMAP> & EMAP,
    std::vector<std::vector<uE2Etype> > & uE2E,
    std::vector< sharptype > & sharp);
  template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedSE>
  IGL_INLINE void sharp_edges(
    const Eigen::MatrixBase<DerivedV> & V,
    const Eigen::MatrixBase<DerivedF> & F,
    const typename DerivedV::Scalar angle,
    Eigen::PlainObjectBase<DerivedSE> & SE
    );
}

#ifndef IGL_STATIC_LIBRARY
#  include "sharp_edges.cpp"
#endif

#endif 

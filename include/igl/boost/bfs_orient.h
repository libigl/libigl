#ifndef IGL_BFS_ORIENT_H
#define IGL_BFS_ORIENT_H
#include <Eigen/Core>
#include <igl/igl_inline.h>

namespace igl
{
  // Consistently orient faces in orientable patches using BFS
  //
  // F = bfs_orient(F,V);
  //
  // Inputs:
  //  F  #F by 3 list of faces
  // Outputs:
  //  FF  #F by 3 list of faces (OK if same as F)
  //  C  #F list of component ids
  //
  //
  template <typename DerivedF, typename DerivedFF, typename DerivedC>
  void bfs_orient(
    const Eigen::PlainObjectBase<DerivedF> & F,
    Eigen::PlainObjectBase<DerivedFF> & FF,
    Eigen::PlainObjectBase<DerivedC> & C);
};
#ifdef IGL_HEADER_ONLY
#  include "bfs_orient.cpp"
#endif

#endif

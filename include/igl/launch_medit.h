#ifndef IGL_LAUNCH_MEDIT_H
#define IGL_LAUNCH_MEDIT_H
#include "igl_inline.h"

#include <Eigen/Dense>

namespace igl 
{
  // Writes the tetmesh in (V,T,F) to a temporary file, opens it with medit
  // (forking with a system call) and returns
  //
  //
  // Templates:
  //   DerivedV  real-value: i.e. from MatrixXd
  //   DerivedT  integer-value: i.e. from MatrixXi
  //   DerivedF  integer-value: i.e. from MatrixXi
  // Inputs:
  //   V  double matrix of vertex positions  #V by 3
  //   T  #T list of tet indices into vertex positions
  //   F  #F list of face indices into vertex positions
  // Returns returned value of system call (probably not useful because of the
  //   fork)
  template <typename DerivedV, typename DerivedT, typename DerivedF>
  IGL_INLINE int launch_medit(
    const Eigen::MatrixBase<DerivedV> & V, 
    const Eigen::MatrixBase<DerivedT> & T,
    const Eigen::MatrixBase<DerivedF> & F);
}

#ifdef IGL_HEADER_ONLY
#  include "launch_medit.cpp"
#endif

#endif


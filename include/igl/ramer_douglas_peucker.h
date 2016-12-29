#ifndef IGL_RAMER_DOUGLAS_PEUCKER_H
#define IGL_RAMER_DOUGLAS_PEUCKER_H
#include "igl_inline.h"
#include <Eigen/Core>
namespace igl
{
  // Ramer-Douglas-Peucker piecewise-linear curve simplification.
  //
  // Inputs:
  //   P  #P by dim ordered list of vertices along the curve
  //   tol  tolerance (maximal euclidean distance allowed between the new line
  //     and a vertex)
  // Outputs:
  //   S  #S by dim ordered list of points along the curve
  //   J  #S list of indices into P so that S = P(J,:)
  template <typename DerivedP, typename DerivedS, typename DerivedJ>
  IGL_INLINE void ramer_douglas_peucker(
    const Eigen::MatrixBase<DerivedP> & P,
    const typename DerivedP::Scalar tol,
    Eigen::PlainObjectBase<DerivedS> & S,
    Eigen::PlainObjectBase<DerivedJ> & J);
}
#ifndef IGL_STATIC_LIBRARY
#  include "ramer_douglas_peucker.cpp"
#endif
#endif

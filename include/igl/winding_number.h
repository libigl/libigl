#ifndef IGL_WINDING_NUMBER_H
#define IGL_WINDING_NUMBER_H
#include "igl_inline.h"
#include <Eigen/Core>

// Minimum number of iterms per openmp thread
#ifndef IGL_WINDING_NUMBER_OMP_MIN_VALUE
#  define IGL_WINDING_NUMBER_OMP_MIN_VALUE 1000
#endif
namespace igl
{
  // WINDING_NUMBER Compute the sum of solid angles of a triangle/tetrahedron
  // described by points (vectors) V
  //
  // Templates:
  //   dim  dimension of input
  // Inputs:
  //  V  n by 3 list of vertex positions
  //  n  number of mesh vertices
  //  F  #F by 3 list of triangle indices, minimum index is 0
  //  m  number of faces
  //  O  no by 3 list of origin positions
  //  no  number of origins
  // Outputs:
  //  S  no by 1 list of winding numbers
  //
  // 3d
  IGL_INLINE void winding_number(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & O,
    Eigen::VectorXd & W);
  template <typename DerivedF>
  IGL_INLINE void winding_number_3(
    const double * V,
    const int n,
    const DerivedF * F,
    const int m,
    const double * O,
    const int no,
    double * S);
  //// Only one evaluation origin
  //template <typename DerivedF>
  //IGL_INLINE void winding_number_3(
  //  const double * V,
  //  const int n,
  //  const DerivedF * F,
  //  const int m,
  //  const double * O,
  //  double * S);
  // 2d
  template <typename DerivedF>
  IGL_INLINE void winding_number_2(
    const double * V,
    const int n,
    const DerivedF * F,
    const int m,
    const double * O,
    const int no,
    double * S);
}

#ifndef IGL_STATIC_LIBRARY
#  include "winding_number.cpp"
#endif

#endif

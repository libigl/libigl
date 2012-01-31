#ifndef IGL_MODE_H
#define IGL_MODE_H
#include "igl_inline.h"
#include <Eigen/Dense>
namespace igl
{
  // Takes mode of coefficients in a matrix along a given dension
  //
  // Templates:
  //   T  should be a eigen matrix primitive type like int or double
  // Inputs:
  //   X  m by n original matrix
  //   d  dension along which to take mode, m or n
  // Outputs:
  //   M  vector containing mode along dension d, if d==1 then this will be a
  //     n-long vector if d==2 then this will be a m-long vector
  template <typename T>
  IGL_INLINE void mode(
    const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & X,
    const int d, 
    Eigen::Matrix<T,Eigen::Dynamic,1> & M);
}

#ifdef IGL_HEADER_ONLY
#  include "mode.cpp"
#endif

#endif

#ifndef IGL_INTERSECT_H
#define IGL_INTERSECT_H
#include <igl/igl_inline.h>

#include <Eigen/Dense>

#ifdef MEX
#  include <mex.h>
#  include <cassert>
#  undef assert
#  define assert( isOK ) ( (isOK) ? (void)0 : (void) mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )
#endif

namespace igl
{
  IGL_INLINE void intersect(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & U,
    const Eigen::MatrixXi & G,
    const bool first_only,
    Eigen::MatrixXi & IF);
}

#ifdef IGL_HEADER_ONLY
#  include "intersect.cpp"
#endif
  
#endif


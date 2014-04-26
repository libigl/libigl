#ifndef IGL_SELFINTERSECT_H
#define IGL_SELFINTERSECT_H
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
  // Optional Parameters
  //   DetectOnly  Only compute IF, leave VV and FF alone
  struct SelfintersectParam
  {
    bool detect_only;
    bool first_only;
    SelfintersectParam():detect_only(false),first_only(false){};
  };
  
  // Given a triangle mesh (V,F) compute a new mesh (VV,FF) which is the same as
  // (V,F) except that any self-intersecting triangles in (V,F) have been
  // subdivided (new vertices and face created) so that the self-intersection
  // contour lies exactly on edges in (VV,FF). New vertices will appear in
  // original faces or on original edges. New vertices on edges are "merged" only
  // across original faces sharing that edge. This means that if the input
  // triangle mesh is a closed manifold the output will be too.
  //
  // Inputs:
  //   V  #V by 3 list of vertex positions
  //   F  #F by 3 list of triangle indices into V
  //   params  struct of optional parameters
  // Outputs:
  //   VV  #VV by 3 list of vertex positions
  //   FF  #FF by 3 list of triangle indices into V
  //   IF  #intersecting face pairs by 2  list of intersecting face pairs,
  //     indexing F
  IGL_INLINE void selfintersect(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const SelfintersectParam & params,
    Eigen::MatrixXd & VV,
    Eigen::MatrixXi & FF,
    Eigen::MatrixXi & IF);
}

#ifdef IGL_HEADER_ONLY
#  include "selfintersect.cpp"
#endif
  
#endif

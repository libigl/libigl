#ifndef IGL_TRIANGULATE_H
#define IGL_TRIANGULATE_H
#include "igl_inline.h"

#ifndef IGL_NO_EIGEN
#  include <Eigen/Core>
#endif
#include <vector>

namespace igl 
{
  // Triangulate a general polygonal mesh into a triangle mesh.
  //
  // Inputs:
  //   vF  list of polygon index lists
  // Outputs:
  //   F  eigen int matrix #F by 3
  //
  // Example:
  //   vector<vector<double > > vV;
  //   vector<vector<int > > vF;
  //   read("poly.obj",vV,vF);
  //   MatrixXd V;
  //   MatrixXi F;
  //   list_to_matrix(vV,V);
  //   triangulate(vF,F);
  template <typename Index, typename DerivedF>
  IGL_INLINE void triangulate(
    const std::vector<std::vector<Index> > & vF,
    Eigen::PlainObjectBase<DerivedF>& F);
}

#ifdef IGL_HEADER_ONLY
#  include "triangulate.cpp"
#endif

#endif


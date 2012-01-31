#ifndef IGL_VF_H
#define IGL_VF_H
#include "igl_inline.h"

#include <Eigen/Dense>
#include <vector>

namespace igl 
{
  // Constructs the vertex-face topology of a given mesh (V,F)
  // Inputs:
  //   V  #V by 3 list of vertex coordinates
  //   F  #F by dim list of mesh faces (must be triangles)
  // Outputs: 
  // 
  //
  // See also: edges, cotmatrix, diag, vv
    
  template <typename T, typename S>
  IGL_INLINE void vf( 
    const Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic> & V, 
    const Eigen::MatrixXi & F, 
    std::vector<std::vector<T> >& VF, std::vector<std::vector<T> >& VFi);
}

#ifdef IGL_HEADER_ONLY
#  include "vf.cpp"
#endif

#endif

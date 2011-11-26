#ifndef IGL_VF_H
#define IGL_VF_H

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Dense>
#include <Eigen/Sparse>

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
    
  template <typename T>
  inline void vf( 
    const Eigen::MatrixXd & V, 
    const Eigen::MatrixXi & F, 
    vector<vector<T> >& VF, vector<vector<T> >& VFi);
}

// Implementation
#include "verbose.h"

template <typename T>
inline void igl::vf(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F, 
  vector<vector<T> >& VF, vector<vector<T> >& VFi)
{
  VF.clear();
  VFi.clear();
  
  VF.resize(V.rows());
  VFi.resize(V.rows());
  
  for(int fi=0; fi<F.rows(); ++fi)
  {
    for(int i = 0; i < 3; ++i)
    {
      VF[F(fi,i)].push_back(fi);
      VFi[F(fi,i)].push_back(i);
    }
  }
}

#endif

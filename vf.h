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
  //   A  vector<vector<int>> of face indices, each row i corresponding to V(i,:)
  //
  // See also: edges, cotmatrix, diag, vv
    
  template <typename T>
  inline void vf( 
    const Eigen::MatrixXd & V, 
    const Eigen::MatrixXi & F, 
    vector<vector<T> >& Al);
}

// Implementation
#include "verbose.h"

template <typename T>
inline void igl::vf(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXi & F, 
  vector<vector<T> >& Al)
{
    Al.clear;
    Al.resize(V.rows());
    
  // Loop over faces
  for(int i = 0;i<F.rows();i++)
  {
    // Loop over this face
    for(int j = 0;j<F.cols();j++)
    {
      Al[F(i,j)].push_back();
    }
  }

  A = Eigen::SparseMatrix<T>(dyn_A);
}

#endif

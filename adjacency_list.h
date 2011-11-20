#ifndef IGL_ADJACENCY_MATRIX_H
#define IGL_ADJACENCY_MATRIX_H

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace igl 
{
  // Constructs the graph adjacency list  of a given mesh (V,F)
  // Templates:
  //   T  should be a eigen sparse matrix primitive type like int or double
  // Inputs:
  //   F  #F by dim list of mesh faces (must be triangles)
  // Outputs: 
  //   A  vector<vector<T> > containing at row i the adjacent vertices of vertex i
  //
  // Example:
  //   // Mesh in (V,F)
  //   vector<vector<double> > A;
  //   adjacency_list(F,A);
  //
  // See also: edges, cotmatrix, diag
  template <typename T>
  inline void adjacency_list(
                             const Eigen::MatrixXi & F, 
                             std::vector<std::vector<T> >& A
                             );
}

// Implementation
#include "verbose.h"

template <typename T>
inline void igl::adjacency_list(
  const Eigen::MatrixXi & F, 
  std::vector<std::vector<T> >& A)
{
  A.clear(); 
  A.resize(F.maxCoeff()+1);

  // Loop over faces
  for(int i = 0;i<F.rows();i++)
  {
    // Loop over this face
    for(int j = 0;j<F.cols();j++)
    {
      // Get indices of edge: s --> d
      int s = F(i,j);
      int d = F(i,(j+1)%F.cols());
      A[s].push_back(d);
      A[d].push_back(s);
    }
  }
  
  // Remove duplicates
  for(int i=0; i<A.size();++i)
  {
    std::sort(A[i].begin(), A[i].end());
    A[i].erase(std::unique(A[i].begin(), A[i].end()), A[i].end());
  }
}

#endif

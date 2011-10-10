#ifndef IGL_EDGES_H
#define IGL_EDGES_H

#include <Eigen/Dense>

namespace igl
{
  // Constructs a list of unique edges represented in a given mesh (V,F)
  // Templates:
  //   T  should be a eigen sparse matrix primitive type like int or double
  // Inputs:
  //   F  #F by 3 list of mesh faces (must be triangles)
  //   or
  //   T  #T x 4  matrix of indices of tet corners
  // Outputs:
  //   E #E by 2 list of edges in no particular order
  //
  // See also: adjacency_matrix
  inline void edges(
    const Eigen::MatrixXi& F, 
    Eigen::MatrixXi& E);
}

// Implementation
#include <map>
#include "adjacency_matrix.h"

inline void igl::edges(
  const Eigen::MatrixXi& F, 
  Eigen::MatrixXi& E)
{
  // build adjacency matrix
  Eigen::SparseMatrix<int> A;
  igl::adjacency_matrix(F,A);
  // Number of non zeros should be twice number of edges
  assert(A.nonZeros()%2 == 0);
  // Resize to fit edges
  E.resize(A.nonZeros()/2,2);
  int i = 0;
  // Iterate over outside
  for(int k=0; k<A.outerSize(); ++k)
  {
    // Iterate over inside
    for(Eigen::SparseMatrix<int>::InnerIterator it (A,k); it; ++it)
    {
      // only add edge in one direction
      if(it.row()<it.col())
      {
        E(i,0) = it.row();
        E(i,1) = it.col();
        i++;
      }
    }
  }
}

#endif

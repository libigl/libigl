#include "edges.h"

#include "adjacency_matrix.h"

IGL_INLINE void igl::edges( const Eigen::MatrixXi& F, Eigen::MatrixXi& E)
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

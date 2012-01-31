#ifndef IGL_ADJACENCY_LIST_H
#define IGL_ADJACENCY_LIST_H
#include "igl_inline.h"

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Dense>
#include <Eigen/Sparse>
//#include <plot_vector.h>
namespace igl 
{
  // Constructs the graph adjacency list of a given mesh (V,F)
  // Templates:
  //   T  should be a eigen sparse matrix primitive type like int or double
  // Inputs:
  //   F       #F by dim list of mesh faces (must be triangles)
  //   sorted  flag that indicates if the list should be sorted counter-clockwise
  // Outputs: 
  //   A  vector<vector<T> > containing at row i the adjacent vertices of vertex i
  //
  // Example:
  //   // Mesh in (V,F)
  //   vector<vector<double> > A;
  //   adjacency_list(F,A);
  //
  // See also: edges, cotmatrix, diag
  template <typename T, typename M>
  IGL_INLINE void adjacency_list(
    const M & F, 
    std::vector<std::vector<T> >& A,
    bool sorted = false);
}

#ifdef IGL_HEADER_ONLY
#  include "adjacency_list.cpp"
#endif

#endif

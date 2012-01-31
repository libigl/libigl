#ifndef IGL_SORT_H
#define IGL_SORT_H
#include "igl_inline.h"

#include <vector>
#include <Eigen/Core>
namespace igl
{

  // Sort the elements of a matrix X along a given dimension like matlabs sort
  // function
  //
  // Templates:
  //   T  should be a eigen matrix primitive type like int or double
  // Inputs:
  //   X  m by n matrix whose entries are to be sorted
  //   dim  dimensional along which to sort:
  //     1  sort each column (matlab default)
  //     2  sort each row
  //   ascending  sort ascending (true, matlab default) or descending (false)
  // Outputs:
  //   Y  m by n matrix whose entries are sorted
  //   IX  m by n matrix of indices so that if dim = 1, then in matlab notation
  //     for j = 1:n, Y(:,j) = X(I(:,j),j); end
  template <typename T>
  IGL_INLINE void sort(
    const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & X,
    const int dim,
    const bool ascending,
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & Y,
    Eigen::MatrixXi & IX);

  // Act like matlab's [Y,I] = SORT(X) for std library vectors
  // Templates:
  //   T  should be a class that implements the '<' comparator operator
  // Input:
  //   unsorted  unsorted vector
  //   ascending  sort ascending (true, matlab default) or descending (false)
  // Output:
  //   sorted     sorted vector, allowed to be same as unsorted
  //   index_map  an index map such that sorted[i] = unsorted[index_map[i]]
  template <class T>
  IGL_INLINE void sort(
      const std::vector<T> &unsorted,
      const bool ascending,
      std::vector<T> &sorted,
      std::vector<size_t> &index_map);
}

#ifdef IGL_HEADER_ONLY
#  include "sort.cpp"
#endif

#endif

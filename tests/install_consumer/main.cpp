#ifdef _WIN32
// <windows.h> defines min/max macros unless NOMINMAX (exported by libigl) is set.
#include <windows.h>
#endif
#include <igl/cotmatrix.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <algorithm>
#include <iostream>

int main()
{
  Eigen::MatrixXd V(3, 3);
  V << 0, 0, 0,
       1, 0, 0,
       0, 1, 0;
  Eigen::MatrixXi F(1, 3);
  F << 0, 1, 2;
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);
  const auto n = std::max<Eigen::Index>(L.rows(), L.cols());
  std::cout << "cotmatrix: " << n << "x" << n << ", nnz=" << L.nonZeros() << std::endl;
  return (n == 3 && L.nonZeros() == 9) ? 0 : 1;
}

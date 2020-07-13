#include <test_common.h>
#include <igl/slice.h>
#include <igl/slice_sorted.h>
#include <random>

namespace
{
  Eigen::SparseMatrix<double> generate_random_sparse_matrix(int rows, int cols)
  {
    std::mt19937 gen;
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    using T = Eigen::Triplet<double>;
    std::vector<T> tripletList;
    for (int i = 0; i < rows; ++i)
    {
      for (int j = 0; j < cols; ++j)
      {
        auto v_ij = dist(gen);  // generate random number
        if (v_ij < 0.1)
        {
          tripletList.push_back(T(i, j, v_ij));  // if larger than treshold, insert it
        }
      }
    }
    Eigen::SparseMatrix<double> mat(rows, cols);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());  // create the matrix
    return mat;
  }

}  // namespace

TEST_CASE("slice_sorted: correctness", "[igl]")
{
  constexpr int rows = 1e3;
  constexpr int cols = 1e3;

  Eigen::SparseMatrix<double> M = generate_random_sparse_matrix(rows, cols);

  Eigen::Matrix<int, Eigen::Dynamic, 1> R(rows / 2);
  Eigen::Matrix<int, Eigen::Dynamic, 1> C(cols / 2);
  for (int i = 0; i < rows; i += 2) R[i / 2] = i;
  for (int i = 0; i < cols; i += 2) C[i / 2] = i;

  SECTION("correctness")
  {
    // Check for correctness
    Eigen::SparseMatrix<double> A, B;
    igl::slice(M, R, C, A);
    igl::slice_sorted(M, R, C, B);
    REQUIRE((A - B).norm() == 0);
  }
}

TEST_CASE("slice_sorted: benchmark", "[igl]" IGL_DEBUG_OFF)
{
  constexpr int rows = 1e3;
  constexpr int cols = 1e3;

  Eigen::SparseMatrix<double> M = generate_random_sparse_matrix(rows, cols);

  Eigen::Matrix<int, Eigen::Dynamic, 1> R(rows / 2);
  Eigen::Matrix<int, Eigen::Dynamic, 1> C(cols / 2);
  for (int i = 0; i < rows; i += 2) R[i / 2] = i;
  for (int i = 0; i < cols; i += 2) C[i / 2] = i;

  BENCHMARK("igl::slice") {
    Eigen::SparseMatrix<double> A;
    igl::slice(M, R, C, A);
    return A.norm();
  };

  BENCHMARK("igl::slice_sorted") {
    Eigen::SparseMatrix<double> A;
    igl::slice_sorted(M, R, C, A);
    return A.norm();
  };
}


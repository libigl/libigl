#include <test_common.h>
#include <Eigen/Sparse>
#include <igl/repmat.h>


template <typename T, int majorType>
void checkRepmat(
  Eigen::SparseMatrix<T, majorType> &A,
  Eigen::SparseMatrix<T, majorType> &B,
  const int rows,
  const int cols,
  const int r,
  const int c)
{
  for (int i = 0; i < rows; i++)
  {
    for (int j = 0; j < cols; j++)
    {
      for (int ii = 0; ii < r; ii++)
      {
        for (int jj = 0; jj < c; jj++)
        {
          REQUIRE (A.coeff(i, j) == B.coeff(i + ii * rows, j + jj * cols));
        }
      }
    }
  }
}

template <typename T, int majorType>
void testRepmat(
  const int rows,
  const int cols,
  const int r,
  const int c)
{
  Eigen::SparseMatrix<T, majorType> A;
  A = (Eigen::Matrix<T, Eigen::Dynamic,
    Eigen::Dynamic>().setRandom(rows, cols).sparseView());
  Eigen::SparseMatrix<T, majorType> B;
  igl::repmat(A, r, c, B);

  REQUIRE (B.rows() == r * rows);
  REQUIRE (B.cols() == c * cols);
  checkRepmat<T, majorType>(A, B, rows, cols, r, c);
}

TEST_CASE("repmat: sparse rowMajor", "[igl]")
{
  testRepmat<double, Eigen::RowMajor>(4,  5, 2, 4);
  testRepmat<   int, Eigen::RowMajor>(2,  8, 3, 4);
  testRepmat< float, Eigen::RowMajor>(6, 10, 2, 2);
}

TEST_CASE("repmat: sparse colMajor", "[igl]")
{
  testRepmat<double, Eigen::ColMajor>(4,  5, 3, 5);
  testRepmat<   int, Eigen::ColMajor>(2,  8, 3, 5);
  testRepmat< float, Eigen::ColMajor>(6, 10, 3, 5);
}
 

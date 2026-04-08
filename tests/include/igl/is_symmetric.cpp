
#include <test_common.h>
#include <igl/is_symmetric.h>

TEST_CASE("is_symmetric: sparse", "[igl]")
{
  {
    Eigen::MatrixXd M(3,3);
    M<<1,2,3,
       4,5,6,
       7,8,9;
    Eigen::SparseMatrix<double> S = M.sparseView();
    REQUIRE (!igl::is_symmetric(S));
  }
  {
    Eigen::MatrixXd M(3,3);
    M<<1,2,3,
       2,4,5,
       3,5,6;
    Eigen::SparseMatrix<double> S = M.sparseView();
    REQUIRE (igl::is_symmetric(S));
  }
  {
    // zero matrix
    Eigen::SparseMatrix<double> S(4,4);
    REQUIRE (igl::is_symmetric(S));
  }
  {
    // matrix with explicit zero entries (produced by cancellation)
    Eigen::MatrixXd M(3,3);
    M<<1,2,3,
       2,4,5,
       3,5,6;
    Eigen::SparseMatrix<double> S = M.sparseView();
    // subtract to create explicit stored zeros
    Eigen::SparseMatrix<double> S2 = S - S;
    REQUIRE (igl::is_symmetric(S2));
  }
  {
    Eigen::SparseMatrix<double> S(1,1);
    S.insert(0,0) = 7.0;
    REQUIRE (igl::is_symmetric(S));
  }
  {
    Eigen::SparseMatrix<double> S(1,1);
    REQUIRE (igl::is_symmetric(S));
  }
}

TEST_CASE("is_symmetric: dense", "[igl]")
{
  {
    Eigen::MatrixXd M(3,3);
    M<<1,2,3,4,5,6,7,8,9;
    REQUIRE (!igl::is_symmetric(M));
  }
  {
    Eigen::MatrixXd M(3,3);
    M<<1,2,3,
       2,4,5,
       3,5,6;
    REQUIRE (igl::is_symmetric(M));
  }
  {
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(4,4);
    REQUIRE (igl::is_symmetric(M));
  }
  {
    Eigen::MatrixXd M(3,3);
    M<<1,0,3,
       0,0,0,
       3,0,6;
    REQUIRE (igl::is_symmetric(M));
  }
  {
    Eigen::MatrixXd M(1,1);
    M << 5.0;
    REQUIRE (igl::is_symmetric(M));
  }
  {
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(1,1);
    REQUIRE (igl::is_symmetric(M));
  }
}

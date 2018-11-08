
#include <test_common.h>
#include <igl/is_symmetric.h>

TEST_CASE("is_symmetric: sparse", "[igl]")
{
  {
    Eigen::MatrixXd M(3,3);
    M<<1,2,3,4,5,6,7,8,9;
    Eigen::SparseMatrix<double> S = M.sparseView();
    REQUIRE (!igl::is_symmetric(S));
  }
  {
    Eigen::MatrixXd M(3,3);
    M<<1,2,3,2,4,5,3,5,6;
    Eigen::SparseMatrix<double> S = M.sparseView();
    REQUIRE (!igl::is_symmetric(S));
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
    REQUIRE (!igl::is_symmetric(M));
  }
}

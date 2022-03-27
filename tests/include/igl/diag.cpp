#include <test_common.h>
#include <igl/diag.h>

TEST_CASE("diag: dense-vector-to-sparse", "[igl]")
{
  const Eigen::VectorXd v = (Eigen::VectorXd(3)<<1,2,3).finished();
  Eigen::SparseMatrix<double> X;
  igl::diag(v,X);
  const Eigen::MatrixXd X_exact = 
    (Eigen::MatrixXd(3,3)<<1,0,0,0,2,0,0,0,3).finished();
  test_common::assert_eq(Eigen::MatrixXd(X),X_exact);
}

TEST_CASE("diag: sparse-vector-to-sparse", "[igl]")
{
  const Eigen::SparseVector<double> v = (Eigen::VectorXd(3)<<1,0,3).finished().sparseView();
  Eigen::SparseMatrix<double> X;
  igl::diag(v,X);
  const Eigen::MatrixXd X_exact = 
    (Eigen::MatrixXd(3,3)<<1,0,0,0,0,0,0,0,3).finished();
  test_common::assert_eq(Eigen::MatrixXd(X),X_exact);
}

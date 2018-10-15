#include <test_common.h>
#include <igl/slice_into.h>
#include <igl/LinSpaced.h>

TEST(slice_into, dense_identity)
{
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(10,9);
  Eigen::VectorXi I = igl::LinSpaced<Eigen::VectorXi >(A.rows(),0,A.rows()-1);
  Eigen::VectorXi J = igl::LinSpaced<Eigen::VectorXi >(A.cols(),0,A.cols()-1);
  {
    Eigen::MatrixXd B(I.maxCoeff()+1,J.maxCoeff()+1);
    igl::slice_into(A,I,J,B);
    test_common::assert_eq(A,B);
  }
  {
    Eigen::MatrixXd B(I.maxCoeff()+1,A.cols());
    igl::slice_into(A,I,1,B);
    test_common::assert_eq(A,B);
  }
  {
    Eigen::MatrixXd B(A.rows(),J.maxCoeff()+1);
    igl::slice_into(A,J,2,B);
    test_common::assert_eq(A,B);
  }
}

TEST(slice_into,density_reverse)
{
  {
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(10,9);
    Eigen::VectorXi I = igl::LinSpaced<Eigen::VectorXi >(A.rows(),A.rows()-1,0);
    Eigen::VectorXi J = igl::LinSpaced<Eigen::VectorXi >(A.cols(),0,A.cols()-1);
    Eigen::MatrixXd B(I.maxCoeff()+1,J.maxCoeff()+1);
    igl::slice_into(A,I,J,B);
    // reverse rows (i.e., reverse each column vector)
    Eigen::MatrixXd C = A.colwise().reverse().eval();
    test_common::assert_eq(B,C);
  }
  {
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(10,9);
    Eigen::VectorXi I = igl::LinSpaced<Eigen::VectorXi >(A.rows(),0,A.rows()-1);
    Eigen::VectorXi J = igl::LinSpaced<Eigen::VectorXi >(A.cols(),A.cols()-1,0);
    Eigen::MatrixXd B(I.maxCoeff()+1,J.maxCoeff()+1);
    igl::slice_into(A,I,J,B);
    // reverse cols (i.e., reverse each row vector)
    Eigen::MatrixXd C = A.rowwise().reverse().eval();
    test_common::assert_eq(B,C);
  }
}


TEST(slice_into, sparse_identity)
{
  Eigen::SparseMatrix<double> A = Eigen::MatrixXd::Random(10,9).sparseView();
  Eigen::VectorXi I = igl::LinSpaced<Eigen::VectorXi >(A.rows(),0,A.rows()-1);
  Eigen::VectorXi J = igl::LinSpaced<Eigen::VectorXi >(A.cols(),0,A.cols()-1);
  {
    Eigen::SparseMatrix<double> B(I.maxCoeff()+1,J.maxCoeff()+1);
    igl::slice_into(A,I,J,B);
    test_common::assert_eq(A,B);
  }
  {
    Eigen::SparseMatrix<double> B(I.maxCoeff()+1,A.cols());
    igl::slice_into(A,I,1,B);
    test_common::assert_eq(A,B);
  }
  {
    Eigen::SparseMatrix<double> B(A.rows(),J.maxCoeff()+1);
    igl::slice_into(A,J,2,B);
    test_common::assert_eq(A,B);
  }
}

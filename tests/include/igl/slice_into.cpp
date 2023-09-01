#include <test_common.h>
#include <igl/slice_into.h>
#include <igl/LinSpaced.h>
#include <igl/randperm.h>

TEST_CASE("slice_into: eigen-random", "[igl]")
{
  const int m = 100;
  const int n = 100;
  Eigen::MatrixXd X = Eigen::MatrixXd::Random(m,n);
  Eigen::VectorXi I;
  igl::randperm(m,I);
  I = I.head(m/2).eval();
  Eigen::VectorXi J;
  igl::randperm(n,J);
  J = J.head(n/2).eval();
  {
    Eigen::MatrixXd Z = Eigen::MatrixXd::Random(I.size(),J.size());
    Eigen::MatrixXd Yigl = X;
    igl::slice_into(Z,I,J,Yigl);
    Eigen::MatrixXd Yeigen = X;
    Yeigen(I,J) = Z;
    test_common::assert_eq(Yigl,Yeigen);
  }
  {
    Eigen::MatrixXd Z = Eigen::MatrixXd::Random(I.size(),X.cols());
    Eigen::MatrixXd Yigl = X;
    igl::slice_into(Z,I,1,Yigl);
    Eigen::MatrixXd Yeigen = X;
    Yeigen(I,Eigen::all) = Z;
    test_common::assert_eq(Yigl,Yeigen);
  }
  {
    Eigen::MatrixXd Z = Eigen::MatrixXd::Random(X.rows(),J.size());
    Eigen::MatrixXd Yigl = X;
    igl::slice_into(Z,J,2,Yigl);
    Eigen::MatrixXd Yeigen = X;
    Yeigen(Eigen::all,J) = Z;
    test_common::assert_eq(Yigl,Yeigen);
  }
  
}

TEST_CASE("slice_into: dense_identity", "[igl]")
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

TEST_CASE("slice_into: density_reverse", "[igl]")
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


TEST_CASE("slice_into: sparse_identity", "[igl]")
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

#include <iostream>
#include <igl/matlab_format.h>
TEST_CASE("slice_into: every-other", "[igl]")
{
  Eigen::MatrixXd Af(2,2);
  Af<<
    1,0,
    5,6;
  Eigen::SparseMatrix<double> As = Af.sparseView();
  Eigen::MatrixXd Bf(4,4);
  Bf<<
    0,5,0,3,
    0,6,0,4,
    3,0,1,5,
    4,8,0,0;
  Eigen::SparseMatrix<double> Bs = Bf.sparseView();

  Eigen::VectorXi R(2);
  R<<1,3;
  Eigen::VectorXi C(2);
  C<<1,3;
  igl::slice_into(Af,R,C,Bf);
  igl::slice_into(As,R,C,Bs);
  Eigen::MatrixXd res(4,4);
  res<<
    0,5,0,3,
    0,1,0,0,
    3,0,1,5,
    4,5,0,6;
  test_common::assert_eq(Bf,res);
  test_common::assert_eq(Eigen::MatrixXd(Bs),res);
}

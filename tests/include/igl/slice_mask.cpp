#include <test_common.h>
#include <igl/slice_mask.h>
#include <igl/randperm.h>
#include <igl/find.h>

TEST_CASE("slice_mask/find: random", "[igl]")
{
  const int m = 100;
  const int n = 100;
  Eigen::MatrixXd X = Eigen::MatrixXd::Random(m,n);
  Eigen::VectorXi I;
  igl::randperm(m,I);
  Eigen::VectorXi J;
  igl::randperm(n,J);
  Eigen::Array<bool,Eigen::Dynamic,1> M = 
    I.unaryExpr([](const int x) { return x%2 == 0; }).eval();
  Eigen::Array<bool,Eigen::Dynamic,1> N = 
    J.unaryExpr([](const int x) { return x%2 == 0; }).eval();
  {
    Eigen::MatrixXd Yigl;
    igl::slice_mask(X,M,N,Yigl);
    Eigen::MatrixXd Yfind = X(igl::find(M),igl::find(N));
    test_common::assert_eq(Yigl,Yfind);
  }
  {
    Eigen::MatrixXd Yigl;
    igl::slice_mask(X,M,1,Yigl);
    Eigen::MatrixXd Yfind = X(igl::find(M),Eigen::all);
    test_common::assert_eq(Yigl,Yfind);
  }
  {
    Eigen::MatrixXd Yigl;
    igl::slice_mask(X,N,2,Yigl);
    Eigen::MatrixXd Yfind = X(Eigen::all,igl::find(N));
    test_common::assert_eq(Yigl,Yfind);
  }
}


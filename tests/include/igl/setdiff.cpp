#include <test_common.h>
#include <igl/setdiff.h>

TEST_CASE("setdiff: matrix", "[igl]")
{
  // Base cases
  {
    const Eigen::VectorXi A = (Eigen::VectorXi(4)<<1,2,1,3).finished();
    const Eigen::VectorXi B(0,1);
    Eigen::VectorXi C,IA;
    const Eigen::VectorXi cC = (Eigen::VectorXi(3)<<1,2,3).finished();
    const Eigen::VectorXi cIA = (Eigen::VectorXi(3)<<0,1,3).finished();
    igl::setdiff(A,B,C,IA);
    test_common::assert_eq(C,cC);
    test_common::assert_eq(IA,cIA);
  }
  {
    const Eigen::VectorXi A(0,1);
    const Eigen::VectorXi B = (Eigen::VectorXi(4)<<1,2,1,3).finished();
    Eigen::VectorXi C,IA;
    const Eigen::VectorXi cC(0,1);
    const Eigen::VectorXi cIA(0,1);
    igl::setdiff(A,B,C,IA);
    test_common::assert_eq(C,cC);
    test_common::assert_eq(IA,cIA);
  }

  {
    // Monkey test
    Eigen::VectorXi A(12);
    A = (Eigen::VectorXd::Random(A.size(),1).array().abs()*9).cast<int>();
    Eigen::VectorXi B(12);
    B = (Eigen::VectorXd::Random(B.size(),1).array().abs()*9).cast<int>();
    Eigen::VectorXi C,IA;
    igl::setdiff(A,B,C,IA);
    for(int i = 0;i<C.size();i++)
    {
      REQUIRE (A(IA(i)) == C(i));
      for(int j = 0;j<B.size();j++)
      {
        REQUIRE (C(i) != B(j));
      }
    }
  }
}


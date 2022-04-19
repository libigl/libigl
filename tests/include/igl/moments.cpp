#include <test_common.h>
#include <igl/moments.h>

TEST_CASE("moments: tet", "[igl]" )
{
  const Eigen::MatrixXd V = 
    (Eigen::MatrixXd(4,3)<<0,0,0, 1,0,0, 0,1,0, 0,0,1).finished();
  const Eigen::MatrixXi F = 
    (Eigen::MatrixXi(4,3)<< 0,2,1, 0,3,2, 2,3,1, 0,1,3).finished();
  double m0;
  Eigen::Vector3d m1;
  Eigen::Matrix3d m2;
  igl::moments(V,F,m0,m1,m2);
  const double epsilon = 1e-15;

  double gt_m0 = 1.0/6.0;
  Eigen::Vector3d gt_m1(1./24.,1./24.,1./24.);
  Eigen::Matrix3d gt_m2;
  gt_m2 <<
    1./80.,1./480.,1./480.,
    1./480.,1./80.,1./480.,
    1./480.,1./480.,1./80.;
  REQUIRE(m0 < gt_m0+epsilon);
  REQUIRE(m0+epsilon > gt_m0);
  test_common::assert_near(m1,gt_m1,epsilon);
  test_common::assert_near(m2,gt_m2,epsilon);
}

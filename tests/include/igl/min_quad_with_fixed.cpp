#include <test_common.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/EPS.h>

TEST_CASE("min_quad_with_fixed: dense", "[igl]" )
{
  const Eigen::Matrix<double,3,3> H = (Eigen::Matrix<double,3,3>(3,3)<<62,43,76,43,69,62,76,62,101).finished();
  const Eigen::Matrix<double,3,1> f = (Eigen::Matrix<double,3,1>(3,1)<<9,8,5).finished();
  const Eigen::Matrix<double,1,3> A = (Eigen::Matrix<double,1,3>(1,3)<<1,1,1).finished();
  const Eigen::Matrix<double,1,1> b = (Eigen::Matrix<double,1,1>(1,1)<<2).finished();
  const Eigen::Array<bool,3,1> k = (Eigen::Array<bool,3,1>()<<true,false,false).finished();
  const Eigen::Matrix<double,3,1> bc = (Eigen::Matrix<double,3,1>(3,1)<<1,0,0).finished();
  // Windows needs template args spelled out
  const Eigen::Matrix<double,3,1> x = igl::min_quad_with_fixed<double,3,1>(H,f,k,bc,A,b);
  REQUIRE(abs(x(0)- 1.0)<igl::EPS<double>());
  REQUIRE(abs(x(1)- 1.5)<igl::EPS<double>());
  REQUIRE(abs(x(2)- -.5)<igl::EPS<double>());
}

#include <test_common.h>
#include <igl/predicates/ear_clipping.h>
#include <iostream>

TEST_CASE("ear_clipping: simple poly 1", "[igl/predicates]")
{
  // Example1: simple polygon 
  Eigen::MatrixXd polygon(10,2);
  polygon<<2,-3,4,1,5.5,-2,6,2.5,5,1,4,5,3,0,1,1,1,5,0,0;
  Eigen::VectorXi RT,nR,M;
  Eigen::MatrixXi eF;
  Eigen::MatrixXd nP;
  RT.setZero(polygon.rows());
  igl::predicates::ear_clipping(polygon,RT,M,eF,nP);
  REQUIRE(nP.rows() == 0);

}

TEST_CASE("ear_clipping: simple poly 2", "[igl/predicates]")
{
  using namespace std;
  const Eigen::MatrixXd P = (Eigen::MatrixXd(6,2)<<
                        0,-0.212132034355964,
        0.212132034355964,                 0,
        0.106066017177982, 0.106066017177982,
                        0,                 0,
       -0.106066017177982, 0.106066017177982,
       -0.212132034355964,                 0
      ).finished();
  const Eigen::VectorXi RT = Eigen::VectorXi::Zero(P.rows(),1);
  Eigen::VectorXi I;
  Eigen::MatrixXi eF;
  Eigen::MatrixXd nP;
  igl::predicates::ear_clipping(P,RT,I,eF,nP);
  Eigen::MatrixXi ans(4, 3);
  ans << 0,1,2,0,2,3,5,0,3,5,3,4;
  REQUIRE(ans == eF);
  REQUIRE(nP.rows() == 0);

}

TEST_CASE("ear_clipping: simple poly 3", "[igl/predicates]")
{
  using namespace std;
  const Eigen::MatrixXd P = (Eigen::MatrixXd(3,2)<<
      0, 0, 
      1, 0,
      0, 1
      ).finished();
  const Eigen::VectorXi RT = Eigen::VectorXi::Zero(P.rows(),1);
  Eigen::VectorXi I;
  Eigen::MatrixXi eF;
  Eigen::MatrixXd nP;
  igl::predicates::ear_clipping(P,RT,I,eF,nP);
  Eigen::MatrixXi ans(1, 3);
  ans << 2, 0, 1;
  REQUIRE(ans == eF);
  REQUIRE(nP.rows() == 0);

}

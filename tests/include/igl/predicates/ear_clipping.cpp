#include <test_common.h>
#include <igl/predicates/ear_clipping.h>

TEST_CASE("ear_clipping: simple poly 1", "[igl/predicates]")
{
  // Example1: simple polygon 
  Eigen::MatrixXd polygon(10,2);
  polygon<<2,-3,4,1,5.5,-2,6,2.5,5,1,4,5,3,0,1,1,1,5,0,0;
  Eigen::VectorXi RT,nR,I;
  Eigen::MatrixXi eF;
  RT.setZero(polygon.rows());

  igl::predicates::ear_clipping(polygon,RT,eF,I);
  REQUIRE(I.size() == 0);
  Eigen::MatrixXi eF1;
  REQUIRE(igl::predicates::ear_clipping(polygon,eF1));
  REQUIRE(eF1.rows() == polygon.rows()-2);

  // Check that reverse also works
  Eigen::MatrixXd polygon2 = polygon.colwise().reverse();
  Eigen::MatrixXi eF2;
  REQUIRE(igl::predicates::ear_clipping(polygon2,eF2));
  REQUIRE(eF2.rows() == polygon2.rows()-2);
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
  igl::predicates::ear_clipping(P,RT,eF,I);
  Eigen::MatrixXi ans(4, 3);
  ans << 0,1,2,0,2,3,5,0,3,5,3,4;
  REQUIRE(ans == eF);
  REQUIRE(I.size() ==0);

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
  igl::predicates::ear_clipping(P,RT,eF,I);
  Eigen::MatrixXi ans(1, 3);
  ans << 2, 0, 1;
  REQUIRE(ans == eF);
  REQUIRE(I.size() == 0);
}

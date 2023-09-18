#include <test_common.h>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>

TEST_CASE("writeOFF: quads", "[igl]")
{
  // Cube
  Eigen::MatrixXd V(8,3);
  V <<
    0,0,0,
    1,0,0,
    1,1,0,
    0,1,0,
    0,0,1,
    1,0,1,
    1,1,1,
    0,1,1;
  Eigen::MatrixXi Q(6,4);
  Q <<
    0,1,2,3,
    1,5,6,2,
    5,4,7,6,
    4,0,3,7,
    3,2,6,7,
    1,0,4,5;
  igl::writeOFF("cube.off",V,Q);
  Eigen::MatrixXd rV;
  Eigen::MatrixXi rQ;
  igl::readOFF("cube.off",rV,rQ);
  test_common::assert_eq(V,rV);
  test_common::assert_eq(Q,rQ);
}

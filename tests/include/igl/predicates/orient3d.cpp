#include <test_common.h>
#include <igl/predicates/orient2d.h>
#include <igl/predicates/orient3d.h>
#include <igl/predicates/incircle.h>
#include <igl/predicates/insphere.h>
#include <igl/predicates/exactinit.h>
#include <limits>

// Didn't have the stamina to break the tests into separate files but they
// should be
TEST_CASE("orient3d", "[igl][predicates]") 
{
  using namespace igl::predicates;
  {
    Eigen::MatrixXd A(3,3);
    Eigen::MatrixXd B(3,3);
    Eigen::MatrixXd C(3,3);
    Eigen::MatrixXd D(3,3);
    A<<
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.0;
    B<<
      1.0, 0.0, 0.0,
      1.0, 0.0, 0.0,
      1.0, 0.0, 0.0;
    C<<
      0.0, 1.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 1.0, 0.0;
    D<<
      0.0, 0.0, 1.0,
      0.0, 0.0, -1.0,
      0.0, 0.0, 0.0;
    Eigen::VectorXi R;
    igl::predicates::orient3d(A,B,C,D,R);
    REQUIRE(R.size() == 3);
    REQUIRE(R(0) == int(igl::predicates::Orientation::NEGATIVE));
    REQUIRE(R(1) == int(igl::predicates::Orientation::POSITIVE));
    REQUIRE(R(2) == int(igl::predicates::Orientation::COPLANAR));
  }
  {
    Eigen::MatrixXd A(1,3);
    Eigen::MatrixXd B(1,3);
    Eigen::MatrixXd C(1,3);
    Eigen::MatrixXd D(3,3);
    A<<
      0.0, 0.0, 0.0,
    B<<
      1.0, 0.0, 0.0,
    C<<
      0.0, 1.0, 0.0,
    D<<
      0.0, 0.0, 1.0,
      0.0, 0.0, -1.0,
      0.0, 0.0, 0.0;
    Eigen::VectorXi R;
    igl::predicates::orient3d(A,B,C,D,R);
    REQUIRE(R.size() == 3);
    REQUIRE(R(0) == int(igl::predicates::Orientation::NEGATIVE));
    REQUIRE(R(1) == int(igl::predicates::Orientation::POSITIVE));
    REQUIRE(R(2) == int(igl::predicates::Orientation::COPLANAR));
  }
}


#include <test_common.h>
#include <igl/circumradius.h>

TEST_CASE("circumradius: equilateral-triangle", "[igl]" )
{
  // Equilateral triangle
  Eigen::MatrixXd V(3,2);
  V << 0,0, 1,0, 0.5,sqrt(3.0)/2.0;
  Eigen::MatrixXi F(1,3);
  F << 0,1,2;
  Eigen::VectorXd R;
  igl::circumradius(V,F,R);
  REQUIRE (R.size() == 1);
  REQUIRE (R(0) == Approx(sqrt(3.0)/3.0).margin(1e-15));
  Eigen::MatrixXd C;
  Eigen::MatrixXd B;
  igl::circumradius(V,F,R,C,B);
  REQUIRE (R.size() == 1);
  REQUIRE (R(0) == Approx(sqrt(3.0)/3.0).margin(1e-15));
  REQUIRE (C.rows() == 1);
  REQUIRE (C.cols() == 2);
  REQUIRE (C(0,0) == Approx(0.5).margin(1e-15));
  REQUIRE (C(0,1) == Approx(sqrt(3.0)/6.0).margin(1e-15));
  REQUIRE (B.rows() == 1);
  REQUIRE (B.cols() == 3);
  REQUIRE (B(0,0) == Approx(1.0/3.0).margin(1e-15));
  REQUIRE (B(0,1) == Approx(1.0/3.0).margin(1e-15));
  REQUIRE (B(0,2) == Approx(1.0/3.0).margin(1e-15));
}

TEST_CASE("circumradius: right-triangle", "[igl]" )
{
  // Right triangle
  Eigen::MatrixXd V(3,2);
  V << 0,0, 1,0, 0,1;
  Eigen::MatrixXi F(1,3);
  F << 0,1,2;
  Eigen::VectorXd R;
  igl::circumradius(V,F,R);
  REQUIRE (R.size() == 1);
  REQUIRE (R(0) == Approx(sqrt(2.0)/2.0).margin(1e-15));
  Eigen::MatrixXd C;
  Eigen::MatrixXd B;
  igl::circumradius(V,F,R,C,B);
  REQUIRE (R.size() == 1);
  REQUIRE (R(0) == Approx(sqrt(2.0)/2.0).margin(1e-15));
  REQUIRE (C.rows() == 1);
  REQUIRE (C.cols() == 2);
  REQUIRE (C(0,0) == Approx(0.5).margin(1e-15));
  REQUIRE (C(0,1) == Approx(0.5).margin(1e-15));
  REQUIRE (B.rows() == 1);
  REQUIRE (B.cols() == 3);
  REQUIRE (B(0,0) == Approx(0.0).margin(1e-15));
  REQUIRE (B(0,1) == Approx(0.5).margin(1e-15));
  REQUIRE (B(0,2) == Approx(0.5).margin(1e-15));
}

// Obtuse
TEST_CASE("circumradius: obtuse-triangle", "[igl]" )
{
  Eigen::MatrixXd V(3,2);
  V << 0,0, 4,0, 2,1;
  Eigen::MatrixXi F(1,3);
  F << 0,1,2;
  Eigen::VectorXd R;
  igl::circumradius(V,F,R);
  REQUIRE (R.size() == 1);
  REQUIRE (R(0) == Approx(2.5).margin(1e-15));

  Eigen::MatrixXd C;
  Eigen::MatrixXd B;
  igl::circumradius(V,F,R,C,B);
  REQUIRE (R.size() == 1);
  REQUIRE (R(0) == Approx(2.5).margin(1e-15));
  REQUIRE (C.rows() == 1);
  REQUIRE (C.cols() == 2);
  REQUIRE (C(0,0) == Approx(2).margin(1e-15));
  REQUIRE (C(0,1) == Approx(-1.5).margin(1e-15));
  REQUIRE (B.rows() == 1);
  REQUIRE (B.cols() == 3);
  REQUIRE (B(0,0) == Approx(1.25).margin(1e-15));
  REQUIRE (B(0,1) == Approx(1.25).margin(1e-15));
  REQUIRE (B(0,2) == Approx(-1.5).margin(1e-15));
}

// Tetrahedra test cases

TEST_CASE("circumradius: equilateral-tetrahedra", "[igl]" )
{
  // Tetrahedra
  Eigen::MatrixXd V(4,3);
  V << 0,0,0, 1,0,0, 0.5,sqrt(3.0)/2.0,0, 0.5,sqrt(3.0)/6.0,sqrt(2.0/3.0);
  Eigen::MatrixXi T(1,4);
  T << 0,1,2,3;
  Eigen::VectorXd R;
  Eigen::MatrixXd C;
  Eigen::MatrixXd B;
  igl::circumradius(V,T,R,C,B);
  REQUIRE (R.size() == 1);
  REQUIRE (R(0) == Approx(sqrt(3.0/8.0)).margin(1e-15));
  REQUIRE (C.rows() == 1);
  REQUIRE (C.cols() == 3);
  REQUIRE (C(0,0) == Approx(0.5).margin(1e-15));
  REQUIRE (C(0,1) == Approx(sqrt(3.0)/6.0).margin(1e-15));
  REQUIRE (C(0,2) == Approx(sqrt(1.0/24.0)).margin(1e-15));
  REQUIRE (B.rows() == 1);
  REQUIRE (B.cols() == 4);
  REQUIRE (B(0,0) == Approx(0.25).margin(1e-15));
  REQUIRE (B(0,1) == Approx(0.25).margin(1e-15));
  REQUIRE (B(0,2) == Approx(0.25).margin(1e-15));
  REQUIRE (B(0,3) == Approx(0.25).margin(1e-15));
}

TEST_CASE("circumradius: right-tetrahedra", "[igl]" )
{
  // Tetrahedra
  Eigen::MatrixXd V(4,3);
  V << 0,0,0, 1,0,0, 0,1,0, 0,0,1;
  Eigen::MatrixXi T(1,4);
  T << 0,1,2,3;
  Eigen::VectorXd R;
  Eigen::MatrixXd C;
  Eigen::MatrixXd B;
  igl::circumradius(V,T,R,C,B);
  REQUIRE (R.size() == 1);
  REQUIRE (R(0) == Approx(sqrt(3.0)/2.0).margin(1e-15));
  REQUIRE (C.rows() == 1);
  REQUIRE (C.cols() == 3);
  REQUIRE (C(0,0) == Approx(0.5).margin(1e-15));
  REQUIRE (C(0,1) == Approx(0.5).margin(1e-15));
  REQUIRE (C(0,2) == Approx(0.5).margin(1e-15));
  REQUIRE (B.rows() == 1);
  REQUIRE (B.cols() == 4);
  REQUIRE (B(0,0) == Approx(-0.5).margin(1e-15));
  REQUIRE (B(0,1) == Approx( 0.5).margin(1e-15));
  REQUIRE (B(0,2) == Approx( 0.5).margin(1e-15));
  REQUIRE (B(0,3) == Approx( 0.5).margin(1e-15));
}

// Obtuse tetrahedron

TEST_CASE("circumradius: obtuse-tetrahedra", "[igl]" )
{
  Eigen::MatrixXd V(4,3);
  V << 0,0,0, 4,0,0, 2,1,0, 2,1,3;
  Eigen::MatrixXi T(1,4);
  T << 0,1,2,3;
  Eigen::VectorXd R;
  Eigen::MatrixXd C;
  Eigen::MatrixXd B;
  igl::circumradius(V,T,R,C,B);
  REQUIRE (R.size() == 1);
  REQUIRE (R(0) == Approx(sqrt(17.0/2.0)).margin(1e-15));
  REQUIRE (C.rows() == 1);
  REQUIRE (C.cols() == 3);
  REQUIRE (C(0,0) == Approx(2).margin(1e-15));
  REQUIRE (C(0,1) == Approx(-1.5).margin(1e-15));
  REQUIRE (C(0,2) == Approx(1.5).margin(1e-15));
  REQUIRE (B.rows() == 1);
  REQUIRE (B.cols() == 4);
  REQUIRE (B(0,0) == Approx(1.25).margin(1e-15));
  REQUIRE (B(0,1) == Approx(1.25).margin(1e-15));
  REQUIRE (B(0,2) == Approx(-2.0).margin(1e-15));
  REQUIRE (B(0,3) == Approx( 0.5).margin(1e-15));
}

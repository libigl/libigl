#include <test_common.h>
#include <igl/euler_characteristic.h>
#include <igl/read_triangle_mesh.h>
#include <igl/triangulated_grid.h>
#include <igl/matlab_format.h>
#include <iostream>

TEST_CASE("euler_characteristic: cube", "[igl]" )
{
  Eigen::MatrixXi F(12,3);
  F <<
    0,1,2,
    0,2,3,
    0,4,5,
    0,5,1,
    0,3,7,
    0,7,4,
    6,5,4,
    6,4,7,
    6,7,3,
    6,3,2,
    6,2,1,
    6,1,5;
  const int chi = igl::euler_characteristic(F);
  REQUIRE(chi == 2);
}

TEST_CASE("euler_characteristic: triangle", "[igl]" )
{
  Eigen::MatrixXi F(1,3);
  F << 0,1,2;
  const int chi = igl::euler_characteristic(F);
  REQUIRE(chi == 1);
}

TEST_CASE("euler_characteristic: square", "[igl]" )
{
  Eigen::MatrixXi F;
  igl::triangulated_grid(3,3,F);
  const int chi = igl::euler_characteristic(F);
  REQUIRE(chi == 1);
}

TEST_CASE("euler_characteristic: torus", "[igl]" )
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(test_common::data_path("TinyTorus.obj"), V, F);
  const int chi = igl::euler_characteristic(F);
  REQUIRE(chi == 0);
}

#include <test_common.h>
#include <igl/voronoi_mass.h>
#include <igl/matlab_format.h>
#include <iostream>

TEST_CASE("voronoi_mass: equilateral-tetrahedra", "[igl]" )
{
  // Tetrahedra
  Eigen::MatrixXd V(4,3);
  V << 0,0,0, 1,0,0, 0.5,sqrt(3.0)/2.0,0, 0.5,sqrt(3.0)/6.0,sqrt(2.0/3.0);
  Eigen::MatrixXi T(1,4);
  T << 0,1,2,3;
  Eigen::VectorXd M;
  igl::voronoi_mass(V,T,M);
  Eigen::VectorXd Mgt(4);
  Mgt <<
  0.0294627825494395,
  0.0294627825494395,
  0.0294627825494395,
  0.0294627825494395;
  test_common::assert_near(M,Mgt,1e-15);
}

TEST_CASE("voronoi_mass: right-tetrahedra", "[igl]" )
{
  // Tetrahedra
  Eigen::MatrixXd V(4,3);
  V << 0,0,0, 1,0,0, 0,1,0, 0,0,1;
  Eigen::MatrixXi T(1,4);
  T << 0,1,2,3;
  Eigen::VectorXd M;
  igl::voronoi_mass(V,T,M);
  Eigen::VectorXd Mgt(4);
  Mgt <<
  0.0833333333333333,
  0.0277777777777778,
  0.0277777777777778,
  0.0277777777777778;
  test_common::assert_near(M,Mgt,1e-15);
}

// Obtuse tetrahedron

TEST_CASE("voronoi_mass: obtuse-tetrahedra", "[igl]" )
{
  Eigen::MatrixXd V(4,3);
  V << 0,0,0, 4,0,0, 2,1,0, 2,1,3;
  Eigen::MatrixXi T(1,4);
  T << 0,1,2,3;
  Eigen::VectorXd M;
  igl::voronoi_mass(V,T,M);
  Eigen::VectorXd Mgt(4);
  Mgt <<
  0.325,
  0.325,
      1,
   0.35;
  test_common::assert_near(M,Mgt,1e-15);
}


#include <test_common.h>
#include <igl/cut_mesh.h>
#include <igl/boundary_loop.h>
#include <igl/vertex_components.h>
#include <igl/edges.h>

TEST_CASE("seperate mesh", "[igl]") {
  
  Eigen::MatrixXd V(9,3);
  V << 0,0,0,
       0,1,0,
       0,2,0,
       1,2,0,
       2,2,0,
       2,1,0,
       2,0,0,
       1,0,0,
       1,1,0;
  Eigen::MatrixXi F(8,3);
  F << 0,1,8,
       1,2,8,
       2,3,8,
       3,4,8,
       4,5,8,
       5,6,8,
       6,7,8,
       7,0,8;
  
  Eigen::MatrixXi C(8,3);
  C << 1,0,0,
       1,0,0,
       1,1,0,
       1,1,1,
       1,1,1,
       1,1,1,
       1,1,1,
       1,0,1;
  C = 1 - C.array();
  Eigen::VectorXi I;
  igl::cut_mesh(V,F,C,I);
  Eigen::VectorXi count;
  igl::vertex_components(F, count);
  REQUIRE(count.maxCoeff() == 2);

}

TEST_CASE("single edge", "[igl]") {
  
  Eigen::MatrixXd V(9,3);
  V << 0,0,0,
       0,1,0,
       0,2,0,
       1,2,0,
       2,2,0,
       2,1,0,
       2,0,0,
       1,0,0,
       1,1,0;
  Eigen::MatrixXi F(8,3);
  F << 0,1,8,
       1,2,8,
       2,3,8,
       3,4,8,
       4,5,8,
       5,6,8,
       6,7,8,
       7,0,8;
  
  Eigen::MatrixXi C(8,3);
  C << 1,0,1,
       1,1,0,
       1,1,1,
       1,1,1,
       1,1,1,
       1,1,1,
       1,1,1,
       1,1,1;
  C = 1 - C.array();
  Eigen::VectorXi I;
  igl::cut_mesh(V,F,C,I);
  Eigen::VectorXi count;
  igl::vertex_components(F, count);
  Eigen::MatrixXi E;
  igl::edges(F, E);
  const auto euler = V.rows() - E.rows() + F.rows();
  CHECK ( 1 == euler );

}
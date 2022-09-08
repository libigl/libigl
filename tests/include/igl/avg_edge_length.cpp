#include <test_common.h>
#include <igl/avg_edge_length.h>
#include <iostream>


TEST_CASE("avg_edge_length: cube", "[igl]")
{
  //The allowed error for this test
  const double epsilon = 1e-15;
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  //This is a cube of dimensions 1.0x1.0x1.0
  igl::read_triangle_mesh(test_common::data_path("cube.obj"), V, F);
  //Create scaled versions of the cube
  double scale = 1.0;
  double huge_scale = 1.0e8;
  double tiny_scale = 1.0e-8;
  Eigen::MatrixXd V_huge = V * huge_scale;
  Eigen::MatrixXd V_tiny = V * tiny_scale;
  //Prepare another mesh with triangles along side diagonals of the cube
  //These triangles are form a regular tetrahedron of side sqrt(2)
  Eigen::MatrixXi F_tet(4,3);
  F_tet << 4,6,1,
            6,4,3,
            4,1,3,
            1,6,3;

  //1. Check avg_edge_length function
  double side_sq = 1.0; //squared lenght of a side
  double diag_sq = 2.0; //squared lenght of a diagonal
  double avg;

  avg = igl::avg_edge_length(V,F);
  REQUIRE (avg == Approx ((12.*sqrt(side_sq) + 6.*sqrt(diag_sq))/(12.+6.)).margin( epsilon));

  //Check the regular tetrahedron
  avg = igl::avg_edge_length(V,F_tet);
  REQUIRE (avg == Approx (sqrt(diag_sq)).margin( epsilon));


  //Scale the cube to have huge sides
  side_sq = huge_scale * huge_scale;  //squared lenght of a side
  diag_sq = 2.0 * side_sq;  //squared lenght of a diagonal
  avg = igl::avg_edge_length(V_huge,F);
  REQUIRE (avg == Approx ((12.*sqrt(side_sq) + 6.*sqrt(diag_sq))/(12.+6.)).margin( epsilon*huge_scale));

  //Check the equilateral triangles
  avg = igl::avg_edge_length(V_huge,F_tet);
  REQUIRE (avg == Approx (sqrt(diag_sq)).margin( epsilon*huge_scale));


  //Scale the cube to have tiny sides
  side_sq = tiny_scale * tiny_scale;  //squared lenght of a side
  diag_sq = 2.0 * side_sq;  //squared lenght of a diagonal
  avg = igl::avg_edge_length(V_tiny,F);
  REQUIRE (avg == Approx ((12.*sqrt(side_sq) + 6.*sqrt(diag_sq))/(12.+6.)).margin( epsilon*tiny_scale));

  //Check the regular tetrahedron
  avg = igl::avg_edge_length(V_tiny,F_tet);
  REQUIRE (avg == Approx (sqrt(diag_sq)).margin( epsilon*tiny_scale));

}

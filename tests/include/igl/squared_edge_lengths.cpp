#include <test_common.h>
#include <igl/squared_edge_lengths.h>
#include <iostream>


TEST_CASE("squared_edge_lengths: cube", "[igl]")
{
  //The allowed error for this test
  const double epsilon = 1e-15;
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  //This is a cube of dimensions 1.0x1.0x1.0
  test_common::load_mesh("cube.obj", V, F);
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

  //1. Check edge_lengths_squared function
  double side_sq = 1.0; //squared lenght of a side
  double diag_sq = 2.0; //squared lenght of a diagonal
  Eigen::MatrixXd L_sq;
  igl::squared_edge_lengths(V,F,L_sq);
  REQUIRE (L_sq.rows() == F.rows());
  REQUIRE (L_sq.cols() == 3);
  //All edges in unit cube measure 1.0 or sqrt(2) in diagonals
  for(int f = 0;f<L_sq.rows();f++)
  {
    //All edge_lengths_squared must be exactly "side_sq" or "diag_sq"
    for(int e = 0;e<3;e++)
        if (L_sq(f,e) > 1.1)
           REQUIRE (L_sq(f,e) == diag_sq);
        else
           REQUIRE (L_sq(f,e) == side_sq);
    //All sides sum exactly side_sq + side_sq + diag_sq
    REQUIRE (side_sq + side_sq + diag_sq == L_sq.row(f).sum());
  }

  //Check the regular tetrahedron
  igl::squared_edge_lengths(V,F_tet,L_sq);
  REQUIRE (L_sq.rows() == F_tet.rows());
  REQUIRE (L_sq.cols() == 3);
  //All edges measure sqrt(2)
  for(int f = 0;f<L_sq.rows();f++)
  {
      //All edge_lengths_squared must be exactly "diag_sq"
    for(int e = 0;e<3;e++)
       REQUIRE (L_sq(f,e) == 2.0);
  }


  //Scale the cube to have huge sides
  side_sq = huge_scale * huge_scale;  //squared lenght of a side
  diag_sq = 2.0 * side_sq;  //squared lenght of a diagonal
  igl::squared_edge_lengths(V_huge,F,L_sq);
  REQUIRE (L_sq.rows() == F.rows());
  REQUIRE (L_sq.cols() == 3);
  for(int f = 0;f<L_sq.rows();f++)
  {
    //All edge_lengths_squared must be exactly side_sq or diag_sq
    for(int e = 0;e<3;e++)
        if (L_sq(f,e) > 1.1*side_sq)
           REQUIRE (L_sq(f,e) == diag_sq);
        else
           REQUIRE (L_sq(f,e) == side_sq);
    //All sides sum exactly side_sq + side_sq + diag_sq
    REQUIRE (side_sq + side_sq + diag_sq == L_sq.row(f).sum());
  }

  //Check the equilateral triangles
  igl::squared_edge_lengths(V_huge,F_tet,L_sq);
  REQUIRE (L_sq.rows() == F_tet.rows());
  REQUIRE (L_sq.cols() == 3);
  //All edges measure sqrt(2)
  for(int f = 0;f<L_sq.rows();f++)
  {
    //All edge_lengths_squared must be exactly "diag_sq"
    for(int e = 0;e<3;e++)
       REQUIRE (L_sq(f,e) == diag_sq);
  }

  //Scale the cube to have tiny sides
  side_sq = tiny_scale * tiny_scale;  //squared lenght of a side
  diag_sq = 2.0 * side_sq;  //squared lenght of a diagonal
  igl::squared_edge_lengths(V_tiny,F,L_sq);
  REQUIRE (L_sq.rows() == F.rows());
  REQUIRE (L_sq.cols() == 3);
  for(int f = 0;f<L_sq.rows();f++)
  {
    //All edge_lengths_squared must be exactly side_sq or diag_sq
    for(int e = 0;e<3;e++)
        if (L_sq(f,e) > 1.1*side_sq)
           REQUIRE (L_sq(f,e) == diag_sq);
        else
           REQUIRE (L_sq(f,e) == side_sq);
    //All sides sum exactly side_sq + side_sq + diag_sq
    REQUIRE (side_sq + side_sq + diag_sq == L_sq.row(f).sum());
  }

  //Check the regular tetrahedron
  igl::squared_edge_lengths(V_tiny,F_tet,L_sq);
  REQUIRE (L_sq.rows() == F_tet.rows());
  REQUIRE (L_sq.cols() == 3);
  //All edges measure sqrt(2)
  for(int f = 0;f<L_sq.rows();f++)
  {
    //All edge_lengths_squared must be exactly "diag_sq"
    for(int e = 0;e<3;e++)
       REQUIRE (L_sq(f,e) == diag_sq);
  }

}

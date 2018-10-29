#include <test_common.h>
#include <igl/edge_lengths.h>
#include <iostream>

TEST_CASE("edge_lengths: cube", "[igl]")
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

  //2. Check edge_lengths
  double side = 1.0;       //lenght of a side
  double diag = sqrt(2.0); //lenght of a diagonal
  Eigen::MatrixXd L;
  igl::edge_lengths(V,F,L);
  REQUIRE (L.rows() == F.rows());
  REQUIRE (L.cols() == 3);
  //All edges in unit cube measure 1.0 or sqrt(2) in diagonals
  for(int f = 0;f<L.rows();f++)
  {
    //All edge_lengths_squared must be exactly "side" or "diag"
    for(int e = 0;e<3;e++)
        if (L(f,e) > 1.1*side)
           REQUIRE (L(f,e) == diag);
        else
           REQUIRE (L(f,e) == side);
    //All sides sum exactly side + side + diag
    REQUIRE (side + side + diag == L.row(f).sum());
  }

  //Check the cube of huge sides
  scale = 1.0e8;
  side = scale;       //lenght of a side
  diag = scale*sqrt(2.0); //lenght of a diagonal
  igl::edge_lengths(V_huge,F,L);
  REQUIRE (L.rows() == F.rows());
  REQUIRE (L.cols() == 3);
  for(int f = 0;f<L.rows();f++)
  {
    //All edge_lengths_squared must be exactly "side" or "diag"
    for(int e = 0;e<3;e++)
        if (L(f,e) > 1.1*side)
           REQUIRE (L(f,e) == diag);
        else
           REQUIRE (L(f,e) == side);
    //All sides sum exactly side + side + diag
    REQUIRE (side + side + diag == Approx (L.row(f).sum()).margin( epsilon));
  }

  //Check the cube of tiny sides
  scale = 1.0e-8;
  side = scale;       //lenght of a side
  diag = scale*sqrt(2.0); //lenght of a diagonal
  igl::edge_lengths(V_tiny,F,L);
  REQUIRE (L.rows() == F.rows());
  REQUIRE (L.cols() == 3);
  for(int f = 0;f<L.rows();f++)
  {
    //All edge_lengths_squared must be exactly "side" or "diag"
    for(int e = 0;e<3;e++)
        if (L(f,e) > 1.1*side)
           REQUIRE (L(f,e) == diag);
        else
           REQUIRE (L(f,e) == side);
    //All sides sum exactly side + side + diag
    REQUIRE (side + side + diag == L.row(f).sum());
  }

}

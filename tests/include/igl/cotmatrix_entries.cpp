#include <test_common.h>
#include <igl/PI.h>
#include <igl/cotmatrix_entries.h>
#include <igl/edge_lengths.h>
#include <igl/EPS.h>

TEST_CASE("cotmatrix_entries: simple", "[igl]")
{
  //The allowed error for this test
  const double epsilon = 1e-15;

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  //This is a cube of dimensions 1.0x1.0x1.0
  igl::read_triangle_mesh(test_common::data_path("cube.obj"), V, F);

  //Prepare another mesh with triangles along side diagonals of the cube
  //These triangles are form a regular tetrahedron of side sqrt(2)
  Eigen::MatrixXi F_tet(4,3);
  F_tet << 4,6,1,
            6,4,3,
            4,1,3,
            1,6,3;

  //1. Check cotmatrix_entries

  Eigen::MatrixXd C1;
  igl::cotmatrix_entries(V,F,C1);
  REQUIRE (C1.rows() == F.rows());
  REQUIRE (C1.cols() == 3);
  //All angles in unit cube measure 45 or 90 degrees
  //Their (half)cotangent must value 0.5 or 0.0
  for(int f = 0;f<C1.rows();f++)
  {
#ifdef IGL_EDGE_LENGTHS_SQUARED_H
      //Hard assert if we have edge_length_squared
      for(int v = 0;v<3;v++)
        if (C1(f,v) > 0.1)
            REQUIRE (C1(f,v) == 0.5);
        else
           REQUIRE (C1(f,v) == 0.0);
       //All cotangents sum 1.0 for those triangles
       REQUIRE (C1.row(f).sum() == 1.0);
#else
      //Soft assert if we have not edge_length_squared
      for(int v = 0;v<3;v++)
        if (C1(f,v) > 0.1)
            REQUIRE (C1(f,v) == Approx (0.5).margin( epsilon));
        else
            REQUIRE (C1(f,v) == Approx (0.0).margin( epsilon));
       //All cotangents sum 1.0 for those triangles
       REQUIRE (C1.row(f).sum() == Approx (1.0).margin( epsilon));
#endif
  }

  //Check the regular tetrahedron
  Eigen::MatrixXd C2;
  igl::cotmatrix_entries(V,F_tet,C2);
  REQUIRE (C2.rows() == F_tet.rows());
  REQUIRE (C2.cols() == 3);
  for(int f = 0;f<C2.rows();f++)
  {
    //Their (half)cotangent must value 0.5 / tan(igl::PI / 3.0)
    for(int v = 0;v<3;v++)
       REQUIRE (C2(f,v) == Approx (0.5 / tan(igl::PI / 3.0)).margin( epsilon));
  }

  //Scale the cube to have huge sides
  Eigen::MatrixXd V_huge = V * 1.0e8;
  igl::cotmatrix_entries(V_huge,F,C1);
  REQUIRE (C1.rows() == F.rows());
  REQUIRE (C1.cols() == 3);
  //All angles still measure 45 or 90 degrees
  //Their (half)cotangent must value 0.5 or 0.0
  for(int f = 0;f<C1.rows();f++)
  {
#ifdef IGL_EDGE_LENGTHS_SQUARED_H
      //Hard assert if we have edge_length_squared
      for(int v = 0;v<3;v++)
        if (C1(f,v) > 0.1)
            REQUIRE (C1(f,v) == 0.5);
        else
           REQUIRE (C1(f,v) == 0.0);
       //All cotangents sum 1.0 for those triangles
       REQUIRE (C1.row(f).sum() == 1.0);
#else
      //Soft assert if we have not edge_length_squared
      for(int v = 0;v<3;v++)
        if (C1(f,v) > 0.1)
            REQUIRE (C1(f,v) == Approx (0.5).margin( epsilon));
        else
            REQUIRE (C1(f,v) == Approx (0.0).margin( epsilon));
       //All cotangents sum 1.0 for those triangles
       REQUIRE (C1.row(f).sum() == Approx (1.0).margin( epsilon));
#endif

  }

  //Check the huge regular tetrahedron
  igl::cotmatrix_entries(V_huge,F_tet,C2);
  REQUIRE (C2.rows() == F_tet.rows());
  REQUIRE (C2.cols() == 3);
  for(int f = 0;f<C2.rows();f++)
  {
    //Their (half)cotangent must value 0.5 / tan(igl::PI / 3.0)
    for(int v = 0;v<3;v++)
       REQUIRE (C2(f,v) == Approx (0.5 / tan(igl::PI / 3.0)).margin( epsilon));
  }

  //Scale the cube to have tiny sides
  Eigen::MatrixXd V_tiny = V * 1.0e-8;
  igl::cotmatrix_entries(V_tiny,F,C1);
  REQUIRE (C1.rows() == F.rows());
  REQUIRE (C1.cols() == 3);
  //All angles still measure 45 or 90 degrees
  //Their (half)cotangent must value 0.5 or 0.0
  for(int f = 0;f<C1.rows();f++)
  {
    for(int v = 0;v<3;v++)
        if (C1(f,v) > 0.1)
           REQUIRE (C1(f,v) == Approx (0.5).margin( epsilon));
        else
           REQUIRE (C1(f,v) == Approx (0.0).margin( epsilon));
    //All cotangents sum 1.0 for those triangles
    REQUIRE (C1.row(f).sum() == Approx (1.0).margin( epsilon));
  }

  //Check the tiny regular tetrahedron
  igl::cotmatrix_entries(V_tiny,F_tet,C2);
  REQUIRE (C2.rows() == F_tet.rows());
  REQUIRE (C2.cols() == 3);
  for(int f = 0;f<C2.rows();f++)
  {
    //Their (half)cotangent must value 0.5 / tan(igl::PI / 3.0)
    for(int v = 0;v<3;v++)
       REQUIRE (C2(f,v) == Approx (0.5 / tan(igl::PI / 3.0)).margin( epsilon));
  }
}// TEST_CASE("cotmatrix_entries: simple", "[igl]")

TEST_CASE("cotmatrix_entries: intrinsic", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  //This is a cube of dimensions 1.0x1.0x1.0
  igl::read_triangle_mesh(test_common::data_path("cube.obj"), V, F);
  Eigen::MatrixXd Cext,Cint;
  // compute C extrinsically
  igl::cotmatrix_entries(V,F,Cext);
  // compute C intrinsically
  Eigen::MatrixXd l;
  igl::edge_lengths(V,F,l);
  igl::cotmatrix_entries(l,Cint);
  test_common::assert_near(Cext,Cint,igl::EPS<double>());
}

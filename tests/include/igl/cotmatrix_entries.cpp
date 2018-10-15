#include <test_common.h>
#include <igl/cotmatrix_entries.h>

TEST(cotmatrix_entries, simple)
{
  //The allowed error for this test
  const double epsilon = 1e-15;

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  //This is a cube of dimensions 1.0x1.0x1.0
  test_common::load_mesh("cube.obj", V, F);

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
  ASSERT_EQ(F.rows(), C1.rows());
  ASSERT_EQ(3, C1.cols());
  //All angles in unit cube measure 45 or 90 degrees
  //Their (half)cotangent must value 0.5 or 0.0
  for(int f = 0;f<C1.rows();f++)
  {
#ifdef IGL_EDGE_LENGTHS_SQUARED_H
      //Hard assert if we have edge_lenght_squared
      for(int v = 0;v<3;v++)
        if (C1(f,v) > 0.1)
            ASSERT_EQ(0.5, C1(f,v));
        else
           ASSERT_EQ(0.0, C1(f,v));
       //All cotangents sum 1.0 for those triangles
       ASSERT_EQ(1.0, C1.row(f).sum());
#else
      //Soft assert if we have not edge_lenght_squared
      for(int v = 0;v<3;v++)
        if (C1(f,v) > 0.1)
            ASSERT_NEAR(0.5, C1(f,v), epsilon);
        else
            ASSERT_NEAR(0.0, C1(f,v), epsilon);
       //All cotangents sum 1.0 for those triangles
       ASSERT_NEAR(1.0, C1.row(f).sum(), epsilon);
#endif
  }

  //Check the regular tetrahedron
  Eigen::MatrixXd C2;
  igl::cotmatrix_entries(V,F_tet,C2);
  ASSERT_EQ(F_tet.rows(), C2.rows());
  ASSERT_EQ(3, C2.cols());
  for(int f = 0;f<C2.rows();f++)
  {
    //Their (half)cotangent must value 0.5 / tan(M_PI / 3.0)
    for(int v = 0;v<3;v++)
       ASSERT_NEAR(0.5 / tan(M_PI / 3.0), C2(f,v), epsilon);
  }

  //Scale the cube to have huge sides
  Eigen::MatrixXd V_huge = V * 1.0e8;
  igl::cotmatrix_entries(V_huge,F,C1);
  ASSERT_EQ(F.rows(), C1.rows());
  ASSERT_EQ(3, C1.cols());
  //All angles still measure 45 or 90 degrees
  //Their (half)cotangent must value 0.5 or 0.0
  for(int f = 0;f<C1.rows();f++)
  {    
#ifdef IGL_EDGE_LENGTHS_SQUARED_H
      //Hard assert if we have edge_lenght_squared
      for(int v = 0;v<3;v++)
        if (C1(f,v) > 0.1)
            ASSERT_EQ(0.5, C1(f,v));
        else
           ASSERT_EQ(0.0, C1(f,v));
       //All cotangents sum 1.0 for those triangles
       ASSERT_EQ(1.0, C1.row(f).sum());
#else
      //Soft assert if we have not edge_lenght_squared
      for(int v = 0;v<3;v++)
        if (C1(f,v) > 0.1)
            ASSERT_NEAR(0.5, C1(f,v), epsilon);
        else
            ASSERT_NEAR(0.0, C1(f,v), epsilon);
       //All cotangents sum 1.0 for those triangles
       ASSERT_NEAR(1.0, C1.row(f).sum(), epsilon);
#endif

  }

  //Check the huge regular tetrahedron
  igl::cotmatrix_entries(V_huge,F_tet,C2);
  ASSERT_EQ(F_tet.rows(), C2.rows());
  ASSERT_EQ(3, C2.cols());
  for(int f = 0;f<C2.rows();f++)
  {
    //Their (half)cotangent must value 0.5 / tan(M_PI / 3.0)
    for(int v = 0;v<3;v++)
       ASSERT_NEAR(0.5 / tan(M_PI / 3.0), C2(f,v), epsilon);
  }

  //Scale the cube to have tiny sides
  Eigen::MatrixXd V_tiny = V * 1.0e-8;
  igl::cotmatrix_entries(V_tiny,F,C1);
  ASSERT_EQ(F.rows(), C1.rows());
  ASSERT_EQ(3, C1.cols());
  //All angles still measure 45 or 90 degrees
  //Their (half)cotangent must value 0.5 or 0.0
  for(int f = 0;f<C1.rows();f++)
  {
    for(int v = 0;v<3;v++)
        if (C1(f,v) > 0.1)
           ASSERT_NEAR(0.5, C1(f,v), epsilon);
        else
           ASSERT_NEAR(0.0, C1(f,v), epsilon);
    //All cotangents sum 1.0 for those triangles
    ASSERT_NEAR(1.0, C1.row(f).sum(), epsilon);
  }

  //Check the tiny regular tetrahedron
  igl::cotmatrix_entries(V_tiny,F_tet,C2);
  ASSERT_EQ(F_tet.rows(), C2.rows());
  ASSERT_EQ(3, C2.cols());
  for(int f = 0;f<C2.rows();f++)
  {
    //Their (half)cotangent must value 0.5 / tan(M_PI / 3.0)
    for(int v = 0;v<3;v++)
       ASSERT_NEAR(0.5 / tan(M_PI / 3.0), C2(f,v), epsilon);
  }

}//TEST(cotmatrix_entries, simple)

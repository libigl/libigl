#include <test_common.h>
#include <igl/cotmatrix.h>

class cotmatrix : public ::testing::TestWithParam<std::string> {};

TEST_P(cotmatrix, constant_in_null_space)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::SparseMatrix<double> L;
  // Load example mesh: GetParam() will be name of mesh file
  test_common::load_mesh(GetParam(), V, F);
  igl::cotmatrix(V,F,L);
  ASSERT_EQ(V.rows(),L.rows());
  ASSERT_EQ(L.rows(),L.cols());
  Eigen::VectorXd C = Eigen::VectorXd::Ones(L.rows());
  Eigen::VectorXd Z = Eigen::VectorXd::Zero(L.rows());
  // ASSERT_EQ(a,b);
  // ASSERT_TRUE(a==b);
  // ASSERT_NEAR(a,b,1e-15)
  ASSERT_LT(((L*C)-(Z)).norm(),1e-12);
}

INSTANTIATE_TEST_CASE_P
(
 all_meshes,
 cotmatrix,
 ::testing::ValuesIn(test_common::all_meshes()),
 test_common::string_test_name
);

TEST(cotmatrix, cube)
{
  //The allowed error for this test
  const double epsilon = 1e-15;

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  //This is a cube of dimensions 1.0x1.0x1.0
  test_common::load_mesh("cube.obj", V, F);

  //Scale the cube to have huge sides
  Eigen::MatrixXd V_huge = V * 1.0e8;

  //Scale the cube to have tiny sides
  Eigen::MatrixXd V_tiny = V * 1.0e-8;

  //Check cotmatrix (Laplacian)
  //The laplacian for the cube is quite singular.
  //Each edge in a diagonal has two opposite angles of 90, with cotangent 0.0 each
  //Each edge in a side has two opposite angle of 45, with (half)cotangen 0.5 each
  //So the cotangent matrix always are (0+0) or (0.5+0.5)
  Eigen::SparseMatrix<double> L1;
  igl::cotmatrix(V,F,L1);
  ASSERT_EQ(V.rows(), L1.rows());
  ASSERT_EQ(V.rows(), L1.cols());
  for(int f = 0;f<L1.rows();f++)
  {
#ifdef IGL_EDGE_LENGTHS_SQUARED_H
    //Hard assert if we have edge_lenght_squared
    ASSERT_EQ(-3.0, L1.coeff(f,f));
    ASSERT_EQ(0.0, L1.row(f).sum());
    ASSERT_EQ(0.0, L1.col(f).sum());
#else
    //Soft assert if we have not edge_lenght_squared
    ASSERT_NEAR(-3.0, L1.coeff(f,f), epsilon);
    ASSERT_NEAR(0.0, L1.row(f).sum(), epsilon);
    ASSERT_NEAR(0.0, L1.col(f).sum(), epsilon);
#endif

  }

  //Same for huge cube.
  igl::cotmatrix(V_huge,F,L1);
  ASSERT_EQ(V.rows(), L1.rows());
  ASSERT_EQ(V.rows(), L1.cols());
  for(int f = 0;f<L1.rows();f++)
  {
    ASSERT_NEAR(-3.0, L1.coeff(f,f), epsilon);
    ASSERT_NEAR(0.0, L1.row(f).sum(), epsilon);
    ASSERT_NEAR(0.0, L1.col(f).sum(), epsilon);
  }

  //Same for tiny cube. we need to use a tolerance this time...
  igl::cotmatrix(V_tiny,F,L1);
  ASSERT_EQ(V.rows(), L1.rows());
  ASSERT_EQ(V.rows(), L1.cols());
  for(int f = 0;f<L1.rows();f++)
  {
    ASSERT_NEAR(-3.0, L1.coeff(f,f), epsilon);
    ASSERT_NEAR(0.0, L1.row(f).sum(), epsilon);
    ASSERT_NEAR(0.0, L1.col(f).sum(), epsilon);
  }

}//TEST(cotmatrix, cube)

TEST(cotmatrix, tetrahedron)
{
  //The allowed error for this test
  const double epsilon = 1e-15;

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  //This is a cube of dimensions 1.0x1.0x1.0
  test_common::load_mesh("cube.obj", V, F);

  //Prepare another mesh with triangles along side diagonals of the cube
  //These triangles are form a regular tetrahedron of side sqrt(2)
  Eigen::MatrixXi F_equi(4,3);
  F_equi << 4,6,1,
            6,4,3,
            4,1,3,
            1,6,3;

  //Scale the cube to have huge sides
  Eigen::MatrixXd V_huge = V * 1.0e8;

  //Scale the cube to have tiny sides
  Eigen::MatrixXd V_tiny = V * 1.0e-8;

  //Check cotmatrix (Laplacian)
  //The laplacian for the cube is quite singular.
  //Each edge in a diagonal has two opposite angles of 90, with cotangent 0.0 each
  //Each edge in a side has two opposite angle of 45, with (half)cotangen 0.5 each
  //So the cotangent matrix always are (0+0) or (0.5+0.5)
  Eigen::SparseMatrix<double> L1;

  //Check the regular tetrahedron of side sqrt(2)
  igl::cotmatrix(V,F_equi,L1);

  ASSERT_EQ(V.rows(), L1.rows());
  ASSERT_EQ(V.rows(), L1.cols());
  for(int f = 0;f<L1.rows();f++)
  {
    //Check the diagonal. Only can value 0.0 for unused vertex or -3 / tan(60)
    if (L1.coeff(f,f) < -0.1)
        ASSERT_NEAR(-3 / tan(M_PI / 3.0), L1.coeff(f,f), epsilon);
    else
        ASSERT_NEAR(0.0, L1.coeff(f,f), epsilon);
#ifdef IGL_EDGE_LENGTHS_SQUARED_H
    //Hard assert if we have edge_lenght_squared
    ASSERT_EQ(0.0, L1.row(f).sum());
    ASSERT_EQ(0.0, L1.col(f).sum());
#else
    //Soft assert if we have not edge_lenght_squared
    ASSERT_NEAR(0.0, L1.row(f).sum(), epsilon);
    ASSERT_NEAR(0.0, L1.col(f).sum(), epsilon);
#endif
  }

  //Check the huge regular tetrahedron
  igl::cotmatrix(V_huge,F_equi,L1);

  ASSERT_EQ(V.rows(), L1.rows());
  ASSERT_EQ(V.rows(), L1.cols());
  for(int f = 0;f<L1.rows();f++)
  {
    //Check the diagonal. Only can value 0.0 for unused vertex or -3 / tan(60)
    if (L1.coeff(f,f) < -0.1)
        ASSERT_NEAR(-3 / tan(M_PI / 3.0), L1.coeff(f,f), epsilon);
    else
        ASSERT_NEAR(0.0, L1.coeff(f,f), epsilon);
    ASSERT_NEAR(0.0, L1.row(f).sum(), epsilon);
    ASSERT_NEAR(0.0, L1.col(f).sum(), epsilon);
  }

  //Check the tiny regular tetrahedron
  igl::cotmatrix(V_tiny,F_equi,L1);

  ASSERT_EQ(V.rows(), L1.rows());
  ASSERT_EQ(V.rows(), L1.cols());
  for(int f = 0;f<L1.rows();f++)
  {
    //Check the diagonal. Only can value 0.0 for unused vertex or -3 / tan(60)
    if (L1.coeff(f,f) < -0.1)
        ASSERT_NEAR(-3 / tan(M_PI / 3.0), L1.coeff(f,f), epsilon);
    else
        ASSERT_NEAR(0.0, L1.coeff(f,f), epsilon);
    ASSERT_NEAR(0.0, L1.row(f).sum(), epsilon);
    ASSERT_NEAR(0.0, L1.col(f).sum(), epsilon);
  }

}//TEST(cotmatrix, tetrahedron)

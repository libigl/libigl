#include <test_common.h>
#include <igl/cotmatrix.h>


TEST_CASE("cotmatrix: constant_in_null_space", "[igl]")
{
  const auto test_case = [](const std::string &param)
  {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::SparseMatrix<double> L;
    // Load example mesh: GetParam() will be name of mesh file
    test_common::load_mesh(param, V, F);
    igl::cotmatrix(V,F,L);
    REQUIRE (L.rows() == V.rows());
    REQUIRE (L.cols() == L.rows());
    Eigen::VectorXd C = Eigen::VectorXd::Ones(L.rows());
    Eigen::VectorXd Z = Eigen::VectorXd::Zero(L.rows());
    // REQUIRE (b == a);
    // REQUIRE (a==b);
    // ASSERT_NEAR(a,b,1e-15)
    REQUIRE (1e-12 > ((L*C)-(Z)).norm());
  };

  test_common::run_test_cases(test_common::all_meshes(), test_case);
}

TEST_CASE("cotmatrix: cube", "[igl]")
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
  REQUIRE (L1.rows() == V.rows());
  REQUIRE (L1.cols() == V.rows());
  for(int f = 0;f<L1.rows();f++)
  {
#ifdef IGL_EDGE_LENGTHS_SQUARED_H
    //Hard assert if we have edge_lenght_squared
    REQUIRE (L1.coeff(f,f) == -3.0);
    REQUIRE (L1.row(f).sum() == 0.0);
    REQUIRE (L1.col(f).sum() == 0.0);
#else
    //Soft assert if we have not edge_lenght_squared
    REQUIRE (L1.coeff(f,f) == Approx (-3.0).margin( epsilon));
    REQUIRE (L1.row(f).sum() == Approx (0.0).margin( epsilon));
    REQUIRE (L1.col(f).sum() == Approx (0.0).margin( epsilon));
#endif

  }

  //Same for huge cube.
  igl::cotmatrix(V_huge,F,L1);
  REQUIRE (L1.rows() == V.rows());
  REQUIRE (L1.cols() == V.rows());
  for(int f = 0;f<L1.rows();f++)
  {
    REQUIRE (L1.coeff(f,f) == Approx (-3.0).margin( epsilon));
    REQUIRE (L1.row(f).sum() == Approx (0.0).margin( epsilon));
    REQUIRE (L1.col(f).sum() == Approx (0.0).margin( epsilon));
  }

  //Same for tiny cube. we need to use a tolerance this time...
  igl::cotmatrix(V_tiny,F,L1);
  REQUIRE (L1.rows() == V.rows());
  REQUIRE (L1.cols() == V.rows());
  for(int f = 0;f<L1.rows();f++)
  {
    REQUIRE (L1.coeff(f,f) == Approx (-3.0).margin( epsilon));
    REQUIRE (L1.row(f).sum() == Approx (0.0).margin( epsilon));
    REQUIRE (L1.col(f).sum() == Approx (0.0).margin( epsilon));
  }
}

TEST_CASE("cotmatrix: tetrahedron", "[igl]")
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

  REQUIRE (L1.rows() == V.rows());
  REQUIRE (L1.cols() == V.rows());
  for(int f = 0;f<L1.rows();f++)
  {
    //Check the diagonal. Only can value 0.0 for unused vertex or -3 / tan(60)
    if (L1.coeff(f,f) < -0.1)
        REQUIRE (L1.coeff(f,f) == Approx (-3 / tan(M_PI / 3.0)).margin( epsilon));
    else
        REQUIRE (L1.coeff(f,f) == Approx (0.0).margin( epsilon));
#ifdef IGL_EDGE_LENGTHS_SQUARED_H
    //Hard assert if we have edge_lenght_squared
    REQUIRE (L1.row(f).sum() == 0.0);
    REQUIRE (L1.col(f).sum() == 0.0);
#else
    //Soft assert if we have not edge_lenght_squared
    REQUIRE (L1.row(f).sum() == Approx (0.0).margin( epsilon));
    REQUIRE (L1.col(f).sum() == Approx (0.0).margin( epsilon));
#endif
  }

  //Check the huge regular tetrahedron
  igl::cotmatrix(V_huge,F_equi,L1);

  REQUIRE (L1.rows() == V.rows());
  REQUIRE (L1.cols() == V.rows());
  for(int f = 0;f<L1.rows();f++)
  {
    //Check the diagonal. Only can value 0.0 for unused vertex or -3 / tan(60)
    if (L1.coeff(f,f) < -0.1)
        REQUIRE (L1.coeff(f,f) == Approx (-3 / tan(M_PI / 3.0)).margin( epsilon));
    else
        REQUIRE (L1.coeff(f,f) == Approx (0.0).margin( epsilon));
    REQUIRE (L1.row(f).sum() == Approx (0.0).margin( epsilon));
    REQUIRE (L1.col(f).sum() == Approx (0.0).margin( epsilon));
  }

  //Check the tiny regular tetrahedron
  igl::cotmatrix(V_tiny,F_equi,L1);

  REQUIRE (L1.rows() == V.rows());
  REQUIRE (L1.cols() == V.rows());
  for(int f = 0;f<L1.rows();f++)
  {
    //Check the diagonal. Only can value 0.0 for unused vertex or -3 / tan(60)
    if (L1.coeff(f,f) < -0.1)
        REQUIRE (L1.coeff(f,f) == Approx (-3 / tan(M_PI / 3.0)).margin( epsilon));
    else
        REQUIRE (L1.coeff(f,f) == Approx (0.0).margin( epsilon));
    REQUIRE (L1.row(f).sum() == Approx (0.0).margin( epsilon));
    REQUIRE (L1.col(f).sum() == Approx (0.0).margin( epsilon));
  }
}

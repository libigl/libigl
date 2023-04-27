#include <test_common.h>
#include <igl/massmatrix.h>
#include <igl/doublearea.h>
#include <igl/edges.h>
#include <igl/vertex_triangle_adjacency.h>

TEST_CASE("massmatrix: full", "[igl]" )
{
  const auto test_case = [](const std::string &param)
  {
    const double epsilon = 1e-15;

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    // Load example mesh: GetParam() will be name of mesh file
    igl::read_triangle_mesh(test_common::data_path(param), V, F);
    Eigen::SparseMatrix<double> M;
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_FULL,M);
    REQUIRE(M.rows() == V.rows());
    REQUIRE(M.cols() == V.rows());
    Eigen::MatrixXd dblA;
    igl::doublearea(V,F,dblA);
    REQUIRE(M.sum() == Approx(dblA.sum()/2.).margin(epsilon));
  };

  test_common::run_test_cases(test_common::all_meshes(), test_case);
}

TEST_CASE("massmatrix: barycentric", "[igl]" )
{
  const auto test_case = [](const std::string &param)
  {
    const double epsilon = 1e-15;

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    // Load example mesh: GetParam() will be name of mesh file
    igl::read_triangle_mesh(test_common::data_path(param), V, F);
    Eigen::SparseMatrix<double> M;
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,M);
    REQUIRE(M.rows() == V.rows());
    REQUIRE(M.cols() == V.rows());
    Eigen::MatrixXd dblA;
    igl::doublearea(V,F,dblA);
    REQUIRE(M.sum() == Approx(dblA.sum()/2.).margin(epsilon));
  };

  test_common::run_test_cases(test_common::all_meshes(), test_case);
}

TEST_CASE("massmatrix: cube_barycentric", "[igl]")
{
  //The allowed error for this test
  const double epsilon = 1e-15;

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  //This is a cube of dimensions 1.0x1.0x1.0
  igl::read_triangle_mesh(test_common::data_path("cube.obj"), V, F);

  Eigen::SparseMatrix<double> M;

  //Check the mass matrix of the cube
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,M);

  REQUIRE(M.rows() == V.rows());
  REQUIRE(M.cols() == V.rows());
  //All triangles' areas are 0.5
  //barycentric area for a vertex is {number of adj triangles} * 0.5 / 3
  Eigen::VectorXi adj(8);
  adj << 6,4,4,4,4,4,6,4;
  for(int f = 0;f<M.rows();f++)
  {
    REQUIRE(M.coeff(f,f) == Approx(0.5/3*adj(f)).margin(epsilon));
  }
}

TEST_CASE("massmatrix: cube_full", "[igl]")
{
  //The allowed error for this test
  const double epsilon = 1e-15;

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  //Check the mass matrix of the cube
  igl::read_triangle_mesh(test_common::data_path("cube.obj"), V, F);

  Eigen::SparseMatrix<double> M;

  //Check the regular tetrahedron of side sqrt(2)
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_FULL,M);

  REQUIRE(M.rows() == V.rows());
  REQUIRE(M.cols() == V.rows());
  Eigen::VectorXi adj(8);
  adj << 6,4,4,4,4,4,6,4;
  //All triangles' areas are 0.5
  //full mass matrix on diagnol is {number of adj triangles} * 0.5 / 6
  for(int f=0;f<M.rows();++f)
  {
    REQUIRE(M.coeff(f,f) == Approx(0.5/6*adj(f)).margin(epsilon));
  }
  for (int i=0;i<F.rows();++i)
  {
    for (int j=0;j<3;++j)
    {
      REQUIRE(M.coeff(F(i,j),F(i,(j+1)%3)) == Approx(0.5/6).margin(epsilon));
      REQUIRE(M.coeff(F(i,(j+1)%3),F(i,j)) == Approx(0.5/6).margin(epsilon));
    }
  }
}

#include <test_common.h>
#include <igl/decimate.h>
#include <igl/sort.h>
#include <igl/sortrows.h>
#include <igl/normalize_row_lengths.h>
#include <igl/slice.h>
#include <igl/matlab_format.h>
#include <iostream>

// class decimate : public ::testing::TestWithParam<std::string> {};

TEST_CASE("decimate: hemisphere", "[igl]")
{
  // Load a hemisphere centered at the origin. For each original vertex compute
  // its "perfect normal" (i.e., its position treated as unit vectors).
  // Decimate the model and using the birth indices of the output vertices grab
  // their original "perfect normals" and compare them to their current
  // positions treated as unit vectors. If vertices have not moved much, then
  // these should be similar (mostly this is checking if the birth indices are
  // sane).
  Eigen::MatrixXd V,U;
  Eigen::MatrixXi F,G;
  Eigen::VectorXi J,I;
  // Load example mesh: GetParam() will be name of mesh file
  test_common::load_mesh("hemisphere.obj", V, F);
  // Perfect normals from positions
  Eigen::MatrixXd NV = V.rowwise().normalized();
  // Remove half of the faces
  igl::decimate(V,F,F.rows()/2,U,G,J,I);
  // Expect that all normals still point in same direction as original
  Eigen::MatrixXd NU = U.rowwise().normalized();
  Eigen::MatrixXd NVI;
  igl::slice(NV,I,1,NVI);
  REQUIRE (NU.rows() == NVI.rows());
  REQUIRE (NU.cols() == NVI.cols());
  // Dot product
  Eigen::VectorXd D = (NU.array()*NVI.array()).rowwise().sum();
  Eigen::VectorXd O = Eigen::VectorXd::Ones(D.rows());
  // 0.2 chosen to succeed on 256 face hemisphere.obj reduced to 128 faces
  test_common::assert_near(D,O,0.02);
}

TEST_CASE("decimate: closed", "[igl]")
{
  const auto test_case = [](const std::string &param)
  {
    Eigen::MatrixXd V,U;
    Eigen::MatrixXi F,G;
    Eigen::VectorXi J;
    // Load example mesh: GetParam() will be name of mesh file
    test_common::load_mesh(param, V, F);
    igl::decimate(V,F,0,U,G,J);
    REQUIRE (4 == U.rows());
    REQUIRE (4 == G.rows());
    {
      Eigen::MatrixXi I;
      igl::sort(Eigen::MatrixXi(G),2,true,G,I);
    }
    {
      Eigen::VectorXi I;
      igl::sortrows(Eigen::MatrixXi(G),true,G,I);
    }
    // Tet with sorted faces
    Eigen::MatrixXi T(4,3);
    T<<
      0,1,2,
      0,1,3,
      0,2,3,
      1,2,3;
    test_common::assert_eq(G,T);
  };

  test_common::run_test_cases(test_common::closed_genus_0_meshes(), test_case);
}
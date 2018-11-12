#include <test_common.h>
#include <igl/cotmatrix_intrinsic.h>
#include <igl/cotmatrix.h>
#include <igl/edge_lengths.h>
#include <igl/matlab_format.h>
#include <igl/EPS.h>

TEST_CASE("cotmatrix_intrinsic: periodic", "[igl]")
{
  const Eigen::MatrixXi F = (Eigen::MatrixXi(18,3)<<
    0,3,1,
    3,4,1,
    1,4,2,
    4,5,2,
    2,5,0,
    5,3,0,
    3,6,4,
    6,7,4,
    4,7,5,
    7,8,5,
    5,8,3,
    8,6,3,
    6,0,7,
    0,1,7,
    7,1,8,
    1,2,8,
    8,2,6,
    2,0,6).finished();
  const Eigen::MatrixXd l = (Eigen::MatrixXd(18,3)<<
    0.47140452079103168,0.33333333333333331,0.33333333333333331,
    0.33333333333333331,0.47140452079103168,0.33333333333333331,
    0.47140452079103168,0.33333333333333331,0.33333333333333331,
    0.33333333333333331,0.47140452079103168,0.33333333333333331,
    0.47140452079103168,0.33333333333333337,0.33333333333333331,
    0.33333333333333331,0.47140452079103168,0.33333333333333337,
    0.47140452079103168,0.33333333333333331,0.33333333333333331,
    0.33333333333333331,0.47140452079103168,0.33333333333333331,
    0.47140452079103168,0.33333333333333331,0.33333333333333331,
    0.33333333333333331,0.47140452079103168,0.33333333333333331,
    0.47140452079103168,0.33333333333333337,0.33333333333333331,
    0.33333333333333331,0.47140452079103168,0.33333333333333337,
    0.47140452079103168,0.33333333333333331,0.33333333333333337,
    0.33333333333333337,0.47140452079103168,0.33333333333333331,
    0.47140452079103168,0.33333333333333331,0.33333333333333337,
    0.33333333333333337,0.47140452079103168,0.33333333333333331,
    0.47140452079103173,0.33333333333333337,0.33333333333333337,
    0.33333333333333337,0.47140452079103173,0.33333333333333337).finished();
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix_intrinsic(l,F,L);
  const Eigen::MatrixXd L_d = L;
  const Eigen::MatrixXd L_gt = (Eigen::MatrixXd(9,9)<<
    -4,1,1,1,0,0,1,0,0,
    1,-4,1,0,1,0,0,1,0,
    1,1,-4,0,0,1,0,0,1,
    1,0,0,-4,1,1,1,0,0,
    0,1,0,1,-4,1,0,1,0,
    0,0,1,1,1,-4,0,0,1,
    1,0,0,1,0,0,-4,1,1,
    0,1,0,0,1,0,1,-4,1,
    0,0,1,0,0,1,1,1,-4).finished();
  test_common::assert_near(L_d,L_gt,igl::EPS<double>());
}

TEST_CASE("cotmatrix_intrinsic: manifold_meshes", "[igl]")
{
  auto test_case = [](const std::string &param)
  {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    test_common::load_mesh(param, V, F);
    Eigen::MatrixXd l;
    igl::edge_lengths(V,F,l);
    Eigen::SparseMatrix<double> L,Li;
    igl::cotmatrix(V,F,L);
    igl::cotmatrix_intrinsic(l,F,Li);
    // Augh, we don't have assert_near for sparse matrices...
    // Instead test that bilinear form is near equal for 2 random vectors
    const Eigen::VectorXd u = Eigen::VectorXd::Random(V.rows(),1);
    const Eigen::VectorXd v = Eigen::VectorXd::Random(V.rows(),1);
    const double uv = u.norm()*v.norm();
    REQUIRE (u.dot(Li*v)/uv == Approx (u.dot(L*v)/uv).margin( igl::EPS<double>()));
  };

  test_common::run_test_cases(test_common::manifold_meshes(), test_case);
}


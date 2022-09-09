#include <test_common.h>
#include <igl/remesh_along_isoline.h>
#include <igl/facet_components.h>
#include <Eigen/Sparse>

TEST_CASE("remesh_along_isoline: triangle_mesh", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  //This is a cube of dimensions 1.0x1.0x1.0
  igl::read_triangle_mesh(test_common::data_path("cube.obj"), V, F);
  const double mean_z = V.col(2).mean();

  Eigen::VectorXi C;
  igl::facet_components(F, C);
  const int fc_count_before = C.maxCoeff();

  Eigen::VectorXd S = V.col(2);
  Eigen::MatrixXd U;
  Eigen::MatrixXi G;
  Eigen::VectorXd SU;
  Eigen::VectorXi J;
  Eigen::SparseMatrix<double> BC;
  Eigen::VectorXi L;
  igl::remesh_along_isoline(V,F,S,mean_z,U,G,SU,J,BC,L);

  igl::facet_components(G, C);
  const int fc_count_after = C.maxCoeff();

  // number of face connected components should not change
  REQUIRE(fc_count_before == fc_count_after);
}

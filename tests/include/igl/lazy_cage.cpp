#include <test_common.h>
#include <igl/lazy_cage.h>
#include <igl/read_triangle_mesh.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/signed_distance.h>
#include <Eigen/Core>
#include <limits>

TEST_CASE("lazy_cage: encloses input", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(test_common::data_path("cube.obj"),V,F);
  const double diag = igl::bounding_box_diagonal(V);

  const int num_faces = 40;
  const int grid_size = 24;
  Eigen::MatrixXd CV;
  Eigen::MatrixXi CF;
  double sigma = -1;
  const bool ok = igl::lazy_cage(
    V,F,num_faces,grid_size,0.5*diag,8,CV,CF,sigma);

  REQUIRE(ok);
  // Target reached exactly (or fewer) and a positive offset was found.
  REQUIRE(CF.rows() > 0);
  REQUIRE(CF.rows() <= num_faces);
  REQUIRE(sigma > 0.0);
  REQUIRE(sigma <= 0.5*diag);

  // The cage must strictly enclose the input: every input vertex is inside the
  // cage (negative signed distance).
  Eigen::VectorXd S; Eigen::VectorXi I; Eigen::MatrixXd C,N;
  igl::signed_distance(
    V,CV,CF,igl::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER,
    -std::numeric_limits<double>::infinity(),
     std::numeric_limits<double>::infinity(),S,I,C,N);
  REQUIRE((S.array() > 1e-7*diag).count() == 0);
}

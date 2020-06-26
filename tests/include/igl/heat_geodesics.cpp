#include <test_common.h>
#include <igl/heat_geodesics.h>
#include <igl/upsample.h>
#include <igl/avg_edge_length.h>

TEST_CASE("heat_geodesic: upsampled cube", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  igl::read_triangle_mesh(test_common::data_path("cube.obj"), V, F);
  igl::upsample(V,F,3);

  int vid = 0;
  igl::HeatGeodesicsData<double> data;
  igl::heat_geodesics_precompute(V,F,data);
  Eigen::VectorXd dist;
  igl::heat_geodesics_solve(data, (Eigen::VectorXi(1,1)<<vid).finished(), dist);
  
  double avg_edge = igl::avg_edge_length(V,F);

  // Check the adjacent corners
  for (int i=0; i<6; i++) {
  REQUIRE((V.row(i)-V.row(0)).norm() == Approx(dist(i)).margin(avg_edge));
  }

}
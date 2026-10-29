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

TEST_CASE("heat_geodesic: semidefinite Poisson system", "[igl]")
{
  Eigen::MatrixXd V(5,3);
  V<<0,0,0,
     1,0,0,
     1,1,0,
     0,1,0,
     0.5,0.5,0;
  Eigen::MatrixXi F(4,3);
  F<<0,1,4,
     1,2,4,
     2,3,4,
     3,0,4;
  for(const int dim : {2,3})
  {
    const Eigen::MatrixXd points = V.leftCols(dim);
    igl::HeatGeodesicsData<double> data;
    REQUIRE(igl::heat_geodesics_precompute(points,F,data));
    Eigen::VectorXi gamma(1);
    gamma<<4;
    Eigen::VectorXd D;
    igl::heat_geodesics_solve(data,gamma,D);
    REQUIRE(D.allFinite());
    REQUIRE(D(4) == Approx(0).margin(1e-12));
    // On each triangle, the normalized heat gradient points from its
    // boundary edge to the center; the discrete Poisson solution is 0.5.
    for(int i=0;i<4;++i)
    {
      REQUIRE(D(i) == Approx(0.5).margin(1e-12));
    }
  }
}

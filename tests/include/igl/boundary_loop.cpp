#include <test_common.h>
#include <igl/boundary_loop.h>
#include <vector>
#include <algorithm>
#include <iostream>

TEST_CASE("boundary_loop: cube", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  //This is a cube of dimensions 1.0x1.0x1.0
  test_common::load_mesh("cube.off", V, F);

  //Compute Boundary Loop
  Eigen::VectorXi boundary;
  igl::boundary_loop(F, boundary);

  //The cube has no boundary
  REQUIRE (boundary.size() == 0);
}

TEST_CASE("boundary_loop: bunny", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  //Load the Stanford bunny
  test_common::load_mesh("bunny.off", V, F);

  //Compute list of ordered boundary loops for a manifold mesh
  std::vector<std::vector<int> >boundaries;
  igl::boundary_loop(F, boundaries);

  //Compare our result with known results taken from meshlab
  REQUIRE (boundaries.size() == 5);

  //Compute min, max and sum of boundaries
  size_t boundaryMin=9999999;
  size_t boundaryMax=0;
  size_t boundarySum=0;
  for(size_t i=0; i<boundaries.size(); i++)
  {
      boundarySum += boundaries[i].size();
      boundaryMax = std::max(boundaryMax, boundaries[i].size());
      boundaryMin = std::min(boundaryMin, boundaries[i].size());
  }

  //Total boundary has 223 vertex
  REQUIRE (boundarySum == 223);
  //Largest loop has 80 vertex
  REQUIRE (boundaryMax == 80);
  //Smallest loop has 22 vertex
  REQUIRE (boundaryMin == 22);
}

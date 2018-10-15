#include <test_common.h>
#include <igl/boundary_loop.h>
#include <vector>
#include <algorithm>
#include <iostream>

TEST(boundary_loop, cube)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  //This is a cube of dimensions 1.0x1.0x1.0
  test_common::load_mesh("cube.off", V, F);

  //Compute Boundary Loop
  Eigen::VectorXi boundary;
  igl::boundary_loop(F, boundary);

  //The cube has no boundary
  ASSERT_EQ(0, boundary.size());
}

TEST(boundary_loop, bunny)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  //Load the Stanford bunny
  test_common::load_mesh("bunny.off", V, F);

  //Compute list of ordered boundary loops for a manifold mesh
  std::vector<std::vector<int> >boundaries;
  igl::boundary_loop(F, boundaries);

  //Compare our result with known results taken from meshlab
  ASSERT_EQ(5, boundaries.size());

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
  ASSERT_EQ(223, boundarySum);
  //Largest loop has 80 vertex
  ASSERT_EQ(80, boundaryMax);
  //Smallest loop has 22 vertex
  ASSERT_EQ(22, boundaryMin);
}

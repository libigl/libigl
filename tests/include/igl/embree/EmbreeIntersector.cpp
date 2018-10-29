#include <test_common.h>
#include <igl/embree/EmbreeIntersector.h>

TEST_CASE("EmbreeIntersector: cube", "[igl/embree]")
{
  //The allowed error for this test
  const double epsilon = 1e-6;

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  // This is a cube of dimensions 1.0x1.0x1.0
  test_common::load_mesh("cube.obj", V, F);

  // Initialize embree
  igl::embree::EmbreeIntersector embree;
  embree.init(V.cast<float>(),F.cast<int>());

  const int expected_id[] = {4,8,5,2,7,0};
  const float expected_u[] = {0.5,0.5,0.5,0.5,0.5,0.5};
  const float expected_v[] = {0.5,0.0,0.0,0.0,0.5,0.0};

  // Shoot ray from inside out
  for (int dim=0; dim<6; ++dim)
  {
    Eigen::Vector3f pos(0,0,0);
    Eigen::Vector3f dir(0,0,0);
    // test each dimension, pos and neg
    dir[dim/2] = dim%2 ? -1 : 1;
    igl::Hit hit;
    bool hitP = embree.intersectRay(pos, dir, hit);
    CHECK(hitP);
    REQUIRE(hit.t == Approx(0.5).margin(epsilon));
    REQUIRE(hit.id == expected_id[dim]);
    REQUIRE(hit.u == Approx(expected_u[dim]).margin(epsilon));
    REQUIRE(hit.v == Approx(expected_v[dim]).margin(epsilon));
  }

  // Shoot ray from outside in
  for (int dim=0; dim<6; ++dim)
  {
    Eigen::Vector3f dir(0,0,0);
    // test each dimension, pos and neg
    dir[dim/2] = dim%2 ? 1 : -1;

    Eigen::Vector3f pos = -dir;

    igl::Hit hit;
    bool hitP = embree.intersectRay(pos, dir, hit);
    CHECK(hitP);
    REQUIRE(hit.t == Approx(0.5).margin(epsilon));
    REQUIRE(hit.id == expected_id[dim]);
    REQUIRE(hit.u == Approx(expected_u[dim]).margin(epsilon));
    REQUIRE(hit.v == Approx(expected_v[dim]).margin(epsilon));
  }

  // Rays that miss
  for (int dim=0; dim<6; ++dim)
  {
    Eigen::Vector3f pos(0,0,0);
    Eigen::Vector3f dir(0,0,0);
    // test each dimension, pos and neg
    dir[dim/2] = dim%2 ? -1 : 1;
    pos[(dim/2+1)%3] = dir[dim/2];

    igl::Hit hit;
    bool hitP = embree.intersectRay(pos, dir, hit);
    CHECK_FALSE(hitP);
  }

  // intersect beam
  {
    Eigen::Vector3f pos(-0.5,-0.5,1);
    Eigen::Vector3f dir(0,0,-1);

    igl::Hit hit;
    bool hitP = embree.intersectBeam(pos, dir, hit);
    CHECK(hitP);
    REQUIRE(hit.t == Approx(0.5).margin(epsilon));
    REQUIRE(hit.id == 7);
    REQUIRE(hit.u == Approx(0).margin(epsilon));
    REQUIRE(hit.v == Approx(1).margin(epsilon));
  }

  {
    Eigen::Vector3f pos(0.5,-1,0.5);
    Eigen::Vector3f dir(0,1,0);

    igl::Hit hit;
    bool hitP = embree.intersectBeam(pos, dir, hit);
    CHECK(hitP);
    REQUIRE(hit.t == Approx(0.5).margin(epsilon));
    REQUIRE(hit.id == 2);
    REQUIRE(hit.u == Approx(0).margin(epsilon));
    REQUIRE(hit.v == Approx(0).margin(epsilon));
  }
}


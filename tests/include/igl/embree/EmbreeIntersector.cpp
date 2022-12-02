#include <test_common.h>
#include <igl/embree/EmbreeIntersector.h>

TEST_CASE("EmbreeIntersector: cube", "[igl/embree]")
{
  //The allowed error for this test
  const double epsilon = 1e-6;

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  // This is a cube of dimensions 1.0x1.0x1.0
  igl::read_triangle_mesh(test_common::data_path("cube.obj"), V, F);

  // Initialize embree
  igl::embree::EmbreeIntersector embree;
  embree.init(V.cast<float>(),F.cast<int>());

  // These are not expected to be exact if the hit is on a vertex of edge.
  const int expected_id[] =  {4,8,5,2,7,0};
  const float expected_u[] = {0.5,0.5,0.5,0.5,0.5,0.5};
  const float expected_v[] = {0.5,0.0,0.0,0.0,0.5,0.0};
  Eigen::MatrixXd hit_P(6,3);
  for (int dim=0; dim<6; ++dim)
  {
    hit_P.row(dim) = 
      V.row(F(expected_id[dim],0))*(1.f - expected_u[dim] - expected_v[dim])+
      V.row(F(expected_id[dim],1))*expected_u[dim] + 
      V.row(F(expected_id[dim],2))*expected_v[dim];
  }
  const auto test_hit = [&](const bool hitP, const igl::Hit& hit, const int dim)
  {
    CHECK(hitP);
    if(hitP)
    {
      const Eigen::RowVectorXd hit_p = 
        V.row(F(hit.id,0))*(1.f - hit.u - hit.v) +
        V.row(F(hit.id,1))*hit.u + 
        V.row(F(hit.id,2))*hit.v;
      // hits will be along diagonal edge so expected_id may be different
      test_common::assert_near(hit_P.row(dim),hit_p,epsilon);
    }
  };


  // Shoot ray from inside out
  for (int dim=0; dim<6; ++dim)
  {
    Eigen::Vector3f pos(0,0,0);
    Eigen::Vector3f dir(0,0,0);
    // test each dimension, pos and neg
    dir[dim/2] = dim%2 ? -1 : 1;
    igl::Hit hit;
    bool hitP = embree.intersectRay(pos, dir, hit);
    test_hit(hitP,hit,dim);
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
    test_hit(hitP,hit,dim);
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
    Eigen::Vector3f pos(1.75,0.25,0);
    Eigen::Vector3f dir(-1,0,0);

    igl::Hit hit;
    bool hitP = embree.intersectBeam(pos, dir, hit);
    CHECK(hitP);
    REQUIRE(hit.t == Approx(1.25).margin(epsilon));
    REQUIRE(hit.id == 4);
    REQUIRE(hit.u == Approx(0.5).margin(epsilon));
    REQUIRE(hit.v == Approx(0.25).margin(epsilon));
  }
}


#include <test_common.h>
#include <igl/ray_mesh_intersect.h>
#include <igl/AABB.h>

TEST_CASE("ray_mesh_intersect: one_triangle", "[igl]")
{
  Eigen::MatrixXd V(3,3);
  V.row(0) << 0.0, 0.0, 0.0;
  V.row(1) << 1.0, 0.0, 0.0;
  V.row(2) << 0.5, 1.0, 0.0;
  
  Eigen::MatrixXi F(1,3);
  F.row(0) << 0,1,2;
  
  Eigen::Vector3f source{0.5, 0.5, -1.0};
  Eigen::Vector3f direction{0.0, 0.0, 1.0};
  
  igl::Hit hit;
  REQUIRE(igl::ray_mesh_intersect(source, direction, V, F, hit) == true);
  REQUIRE(hit.t == Approx(1.0));
  
  std::vector<igl::Hit> hits;
  REQUIRE(igl::ray_mesh_intersect(source, direction, V, F, hits) == true);
  REQUIRE(hits.size() == 1);
  REQUIRE(hits.front().t == Approx(1.0));
}


TEST_CASE("ray_mesh_intersect: corner-case", "[igl]")
{
  //       .      //
  //      /|\     //
  //     / | \    //
  //    /  |  \   //
  //   /   |   \  //
  //  .----x----. //      
  //   \   |   /  //               
  //    \  |  /   //               
  //     \ | /    //               
  //      \|/     //                
  //       .      //

  Eigen::MatrixXf vertices(5, 3);
  vertices <<
      -1., 0., 0.,
      0., 1., 0.,
      1., 0., 0.,
      0., -1., 0.,
      0., 0., 0.;

  Eigen::MatrixXi faces(4, 3);
  faces <<
      0, 1, 4,
      1, 2, 4,
      2, 3, 4,
      3, 0, 4;

  igl::AABB<Eigen::MatrixXf, 3> mesh_bvh;

  mesh_bvh.init(vertices, faces);

  for (float eps: {1e-5f, 0.f})
  {
    Eigen::Vector3f origin(eps, eps, 1.f + eps);
    Eigen::Vector3f direction(0.f, 0.f, -1.f);

    std::vector<igl::Hit> hits, hits_bvh;
    bool is_hit = igl::ray_mesh_intersect(origin, direction, vertices, faces, hits);
    bool is_hit_bvh = mesh_bvh.intersect_ray(vertices, faces, origin, direction, hits_bvh);

    REQUIRE (is_hit);
    REQUIRE (is_hit == is_hit_bvh);
    REQUIRE (hits.size() == hits_bvh.size());
  }
}

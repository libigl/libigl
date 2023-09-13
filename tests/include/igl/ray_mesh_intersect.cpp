#include <test_common.h>
#include <igl/ray_mesh_intersect.h>
#include <igl/ray_box_intersect.h>
#include <igl/AABB.h>


TEST_CASE("ray_mesh_intersect: one_triangle", "[igl]")
{
  IGL_PUSH_FPE;
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

  IGL_POP_FPE;
}

TEST_CASE("ray_mesh_intersect: corner-case", "[igl]")
{
  IGL_PUSH_FPE;
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
  IGL_POP_FPE;
}

TEST_CASE("ray_mesh_intersect: corner-case2", "[igl]")
{
  IGL_PUSH_FPE;

  Eigen::MatrixXf vertices(3, 3);
  vertices <<
      -2.891303300857544, 0.7025225162506104, 1.157850384712219,
      -2.870383024215698, 0.7444183230400085, 1.18663215637207,
      -2.890183448791504, 0.7462523579597473, 1.157822966575623;


  Eigen::MatrixXi faces(1, 3);
  faces <<
      0, 2, 1;

  Eigen::Vector3f origin;
  Eigen::Vector3f direction;

  origin << -5.411622047424316, -0.02165498770773411, 0.7916983366012573;
  direction << 0.9475222229957581, 0.2885690927505493, 0.1375846415758133;

  std::vector<igl::Hit> hits, hits_bvh;
  bool is_hit = igl::ray_mesh_intersect(origin, direction, vertices, faces, hits);
  Eigen::AlignedBox3f box;
  box.extend(vertices.row(0).transpose());
  box.extend(vertices.row(1).transpose());
  box.extend(vertices.row(2).transpose());

  float tmin, tmax;
  bool is_hit_box = igl::ray_box_intersect(origin, direction, box, 0.0f, std::numeric_limits<float>::max(), tmin, tmax);
  REQUIRE (is_hit == is_hit_box);
}

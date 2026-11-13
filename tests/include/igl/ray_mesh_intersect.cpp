#include <test_common.h>
#include <igl/ray_mesh_intersect.h>
#include <igl/ray_box_intersect.h>
#include <igl/AABB.h>
#include <cmath>
#include <limits>


TEST_CASE("ray_mesh_intersect: one_triangle", "[igl]")
{
  IGL_PUSH_FPE;
  Eigen::MatrixXd V(3,3);
  V.row(0) << 0.0, 0.0, 0.0;
  V.row(1) << 1.0, 0.0, 0.0;
  V.row(2) << 0.5, 1.0, 0.0;
  
  Eigen::MatrixXi F(1,3);
  F.row(0) << 0,1,2;
  
  Eigen::Vector3d source{0.5, 0.5, -1.0};
  Eigen::Vector3d direction{0.0, 0.0, 1.0};
  
  igl::Hit<double> hit;
  REQUIRE(igl::ray_mesh_intersect(source, direction, V, F, hit) == true);
  REQUIRE(hit.t == Approx(1.0));
  
  std::vector<igl::Hit<double>> hits;
  REQUIRE(igl::ray_mesh_intersect(source, direction, V, F, hits) == true);
  REQUIRE(hits.size() == 1);
  REQUIRE(hits.front().t == Approx(1.0));

  IGL_POP_FPE;
}

// https://github.com/libigl/libigl/issues/2363
#ifdef _WIN32
#define IGL_MAYFAIL_WIN "[!mayfail]"
#else
#define IGL_MAYFAIL_WIN
#endif
TEST_CASE("ray_mesh_intersect: corner-case", "[igl]" IGL_MAYFAIL_WIN)
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

    std::vector<igl::Hit<float>> hits, hits_bvh;
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

  std::vector<igl::Hit<float>> hits, hits_bvh;
  bool is_hit = igl::ray_mesh_intersect(origin, direction, vertices, faces, hits);
  Eigen::AlignedBox3f box;
  box.extend(vertices.row(0).transpose());
  box.extend(vertices.row(1).transpose());
  box.extend(vertices.row(2).transpose());

  float tmin, tmax;
  bool is_hit_box = igl::ray_box_intersect(origin, direction, box, 0.0f, std::numeric_limits<float>::max(), tmin, tmax);
  REQUIRE (is_hit == is_hit_box);
}

// https://github.com/libigl/libigl/issues/2568
TEMPLATE_TEST_CASE("ray_mesh_intersect: rotated vertex hits", "[igl]",
  float, double)
{
  IGL_PUSH_FPE;
  typedef TestType Scalar;
  typedef Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> Vertices;
  typedef Eigen::Matrix<Scalar,3,1> Vector3;
  Vertices original(5,3);
  original << -1,0,0, 0,1,0, 1,0,0, 0,-1,0, 0,0,0;
  Eigen::MatrixXi F(4,3);
  F << 0,1,4, 1,2,4, 2,3,4, 3,0,4;
  for(int axis = 0; axis < 3; ++axis)
  {
    Vertices V(5,3);
    for(int j = 0; j < 3; ++j)
    {
      V.col((j+axis)%3) = original.col(j);
    }
    igl::AABB<Vertices,3> tree;
    tree.init(V,F);
    for(const Scalar side : {Scalar(-1),Scalar(1)})
    {
      for(const Scalar eps : {Scalar(0),Scalar(1e-5)})
      {
        DYNAMIC_SECTION("axis " << axis << " side " << side << " eps "
          << eps)
        {
          Vector3 origin = Vector3::Constant(eps);
          Vector3 direction = Vector3::Zero();
          origin((2+axis)%3) = side*(Scalar(1)+eps);
          direction((2+axis)%3) = -side;
          std::vector<igl::Hit<Scalar>> direct, accelerated;
          REQUIRE(igl::ray_mesh_intersect(origin,direction,V,F,direct));
          REQUIRE(tree.intersect_ray(V,F,origin,direction,accelerated));
          REQUIRE(direct.size() == accelerated.size());
          if(eps == Scalar(0))
          {
            REQUIRE(accelerated.size() == 4);
          }
        }
      }
    }
  }
  IGL_POP_FPE;
}

TEMPLATE_TEST_CASE("ray_box_intersect: parallel boundary slabs", "[igl]",
  float, double)
{
  IGL_PUSH_FPE;
  typedef TestType Scalar;
  typedef Eigen::Matrix<Scalar,1,3> RowVector3;
  typedef Eigen::Matrix<Scalar,3,1> Vector3;
  const Eigen::AlignedBox<Scalar,3> box(Vector3::Zero(),Vector3::Ones());
  for(int axis = 0; axis < 3; ++axis)
  {
    for(const Scalar side : {Scalar(-1),Scalar(1)})
    {
      for(const Scalar zero : {Scalar(0),-Scalar(0)})
      {
        for(const Scalar position : {Scalar(-0.1),Scalar(0),Scalar(0.5),
          Scalar(1),Scalar(1.1)})
        {
          DYNAMIC_SECTION("axis " << axis << " side " << side << " zero "
            << std::signbit(zero) << " position " << position)
          {
            RowVector3 origin = RowVector3::Constant(position);
            RowVector3 direction = RowVector3::Constant(zero);
            origin(axis) = side > 0 ? Scalar(-1) : Scalar(2);
            direction(axis) = side;
            Scalar tmin, tmax;
            const bool expected = position >= 0 && position <= 1;
            const bool hit = igl::ray_box_intersect(origin,direction,box,
              Scalar(0),Scalar(3),tmin,tmax);
            REQUIRE(hit == expected);
            if(hit)
            {
              REQUIRE(tmin == Approx(1));
              REQUIRE(tmax == Approx(2));
              // A finite segment stopping before the box still misses it.
              REQUIRE_FALSE(igl::ray_box_intersect(origin,direction,box,
                Scalar(0),Scalar(0.5),tmin,tmax));
            }
          }
        }
      }
    }
  }
  IGL_POP_FPE;
}

TEMPLATE_TEST_CASE("ray_box_intersect: negative tangent interval", "[igl]",
  float, double)
{
  IGL_PUSH_FPE;
  typedef TestType Scalar;
  typedef Eigen::Matrix<Scalar,1,3> RowVector3;
  typedef Eigen::Matrix<Scalar,3,1> Vector3;
  const Eigen::AlignedBox<Scalar,3> box(Vector3::Zero(),Vector3::Ones());
  for(int axis = 0; axis < 3; ++axis)
  {
    for(const Scalar zero : {Scalar(0),-Scalar(0)})
    {
      DYNAMIC_SECTION("axis " << axis << " zero " << std::signbit(zero))
      {
        RowVector3 origin, direction;
        origin(axis) = Scalar(1.5);
        origin((axis+1)%3) = Scalar(1);
        origin((axis+2)%3) = Scalar(0.5);
        direction(axis) = Scalar(1);
        direction((axis+1)%3) = zero;
        direction((axis+2)%3) = Scalar(-1);
        Scalar tmin, tmax;
        // The only intersection is the corner reached at t = -0.5.
        REQUIRE(igl::ray_box_intersect(origin,direction,box,
          Scalar(-2),Scalar(0.125),tmin,tmax));
        REQUIRE(tmin == Approx(-0.5));
        REQUIRE(tmax == Approx(-0.5));
        REQUIRE_FALSE(igl::ray_box_intersect(origin,direction,box,
          Scalar(0),Scalar(2),tmin,tmax));
        // A parallel outside slab must remain a miss even far behind us.
        origin(axis) = Scalar(1.125);
        origin((axis+1)%3) = Scalar(3);
        origin((axis+2)%3) = Scalar(0.5);
        direction(axis) = zero;
        direction((axis+1)%3) = std::numeric_limits<Scalar>::min();
        direction((axis+2)%3) = zero;
        REQUIRE_FALSE(igl::ray_box_intersect(origin,direction,box,
          -std::numeric_limits<Scalar>::max(),Scalar(0),tmin,tmax));
      }
    }
  }
  IGL_POP_FPE;
}

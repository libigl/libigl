#include <test_common.h>
#include <igl/exact_ray_mesh_intersect.h>
#include <igl/exact_ray_triangle_intersect.h>
#include <igl/exact_ray_triangle_hit_compare.h>
#include <igl/ray_mesh_intersect.h>
#include <igl/read_triangle_mesh.h>
#include <random>
#include <algorithm>
#include <vector>

TEST_CASE("exact_ray_mesh_intersect: one_triangle", "[igl]")
{
  Eigen::MatrixXd V(3,3);
  V.row(0) << 0.0, 0.0, 0.0;
  V.row(1) << 1.0, 0.0, 0.0;
  V.row(2) << 0.5, 1.0, 0.0;
  Eigen::MatrixXi F(1,3);
  F.row(0) << 0,1,2;

  Eigen::Vector3d source{0.5, 0.4, -1.0};
  Eigen::Vector3d dir{0.0, 0.0, 1.0};

  int fid = -1;
  REQUIRE(igl::exact_ray_mesh_intersect(source, dir, V, F, fid));
  REQUIRE(fid == 0);

  igl::Hit<double> hit;
  REQUIRE(igl::exact_ray_mesh_intersect(source, dir, V, F, hit));
  REQUIRE(hit.id == 0);
  REQUIRE(hit.t == Approx(1.0));

  // Ray pointing away: no hit (t < 0).
  Eigen::Vector3d away{0.0, 0.0, -1.0};
  REQUIRE_FALSE(igl::exact_ray_mesh_intersect(source, away, V, F, fid));

  // Ray missing the triangle.
  Eigen::Vector3d off{5.0, 5.0, -1.0};
  REQUIRE_FALSE(igl::exact_ray_mesh_intersect(off, dir, V, F, fid));
}

TEST_CASE("exact_ray_mesh_intersect: ordering", "[igl]")
{
  // Three parallel triangles at z = 3, 1, 2 (deliberately out of index order).
  Eigen::MatrixXd V(9,3);
  const double zs[3] = {3.0, 1.0, 2.0};
  for(int f = 0; f < 3; f++)
  {
    V.row(3*f+0) << 0.0, 0.0, zs[f];
    V.row(3*f+1) << 2.0, 0.0, zs[f];
    V.row(3*f+2) << 0.0, 2.0, zs[f];
  }
  Eigen::MatrixXi F(3,3);
  F << 0,1,2, 3,4,5, 6,7,8;

  Eigen::Vector3d source{0.3, 0.3, -5.0};
  Eigen::Vector3d dir{0.0, 0.0, 1.0};

  int fid = -1;
  REQUIRE(igl::exact_ray_mesh_intersect(source, dir, V, F, fid));
  REQUIRE(fid == 1); // z = 1 is closest

  std::vector<int> fids;
  REQUIRE(igl::exact_ray_mesh_intersect_all(source, dir, V, F, fids));
  REQUIRE(fids.size() == 3);
  REQUIRE(fids[0] == 1); // z = 1
  REQUIRE(fids[1] == 2); // z = 2
  REQUIRE(fids[2] == 0); // z = 3

  std::vector<igl::Hit<double>> hits;
  REQUIRE(igl::exact_ray_mesh_intersect_all(source, dir, V, F, hits));
  REQUIRE(hits.size() == 3);
  REQUIRE(hits[0].t <= hits[1].t);
  REQUIRE(hits[1].t <= hits[2].t);
}

// The exact leaf boolean and comparison predicates must agree with a brute force
// scan, and the AABB traversal must agree with that scan on random meshes.
TEST_CASE("exact_ray_mesh_intersect: tree_matches_bruteforce", "[igl]")
{
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> U(-1,1);
  for(int t = 0; t < 3000; t++)
  {
    const int nf = 1 + (rng()%16);
    Eigen::MatrixXd V(3*nf,3);
    Eigen::MatrixXi F(nf,3);
    for(int f = 0; f < nf; f++)
    {
      for(int k = 0; k < 3; k++){ V.row(3*f+k) << U(rng),U(rng),U(rng); }
      F.row(f) << 3*f, 3*f+1, 3*f+2;
    }
    Eigen::Vector3d s(2*U(rng), 2*U(rng), 2*U(rng)), d(U(rng),U(rng),U(rng));
    if(d.norm() < 1e-6){ continue; }

    // Fixed-size rows so the leaf predicates are called with concrete types.
    const auto R = [&](int i){ return Eigen::RowVector3d(V.row(i)); };
    const auto tri_before = [&](int a, int b)
    {
      return igl::exact_ray_triangle_hit_compare(s,d,
        R(F(a,0)),R(F(a,1)),R(F(a,2)),
        R(F(b,0)),R(F(b,1)),R(F(b,2))) < 0;
    };
    std::vector<int> bf;
    for(int f = 0; f < nf; f++)
    {
      if(igl::exact_ray_triangle_intersect(
           s,d,R(F(f,0)),R(F(f,1)),R(F(f,2))))
      { bf.push_back(f); }
    }
    std::sort(bf.begin(), bf.end(), tri_before);

    std::vector<int> tr;
    igl::exact_ray_mesh_intersect_all(s,d,V,F,tr);

    // Same set of hit faces.
    std::vector<int> a = bf, b = tr;
    std::sort(a.begin(),a.end()); std::sort(b.begin(),b.end());
    REQUIRE(a == b);

    // First hit consistent (equal exact distance to brute force front).
    int fid = -1;
    const bool h = igl::exact_ray_mesh_intersect(s,d,V,F,fid);
    REQUIRE(h == !bf.empty());
    if(h)
    {
      const int c = igl::exact_ray_triangle_hit_compare(s,d,
        R(F(fid,0)),R(F(fid,1)),R(F(fid,2)),
        R(F(bf.front(),0)),R(F(bf.front(),1)),R(F(bf.front(),2)));
      REQUIRE(c == 0);
    }
  }
}

// On a real, well-separated mesh the exact classification should agree with the
// ordinary double-precision Moller-Trumbore reference.
TEST_CASE("exact_ray_mesh_intersect: agrees_with_double", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(test_common::data_path("cube.obj"), V, F);

  Eigen::RowVector3d bmin = V.colwise().minCoeff();
  Eigen::RowVector3d bmax = V.colwise().maxCoeff();
  Eigen::RowVector3d c = 0.5*(bmin+bmax);
  const double r = (bmax-bmin).norm();

  std::mt19937 rng(11);
  std::uniform_real_distribution<double> U(-1,1);
  int compared = 0;
  for(int i = 0; i < 500; i++)
  {
    Eigen::Vector3d p(U(rng),U(rng),U(rng));
    if(p.norm() < 1e-3){ continue; }
    p.normalize();
    Eigen::Vector3d source = c.transpose() + r*p; // outside, looking in
    Eigen::Vector3d dir = (c.transpose() - source);

    std::vector<igl::Hit<double>> ref;
    igl::ray_mesh_intersect(source, dir, V, F, ref);
    int fid = -1;
    const bool got = igl::exact_ray_mesh_intersect(source, dir, V, F, fid);
    REQUIRE(got == (ref.size() > 0));
    if(got)
    {
      igl::Hit<double> hit;
      igl::exact_ray_mesh_intersect(source, dir, V, F, hit);
      REQUIRE(hit.t == Approx(ref.front().t).margin(1e-9));
      compared++;
    }
  }
  REQUIRE(compared > 0);
}

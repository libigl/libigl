#include <test_common.h>
#include <igl/lazy_cage.h>
#include <igl/read_triangle_mesh.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/signed_distance.h>
#include <igl/centroid.h>
#include <igl/find_intersections.h>
#include <Eigen/Core>
#include <limits>

namespace
{
  // A valid cage must enclose the input (all input vertices inside) and not
  // intersect it.
  void require_valid_cage(
    const Eigen::MatrixXd & V, const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & CV, const Eigen::MatrixXi & CF)
  {
    const double diag = igl::bounding_box_diagonal(V);
    // Encloses input.
    Eigen::VectorXd S; Eigen::VectorXi I; Eigen::MatrixXd C,N;
    igl::signed_distance(
      V,CV,CF,igl::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER,
      -std::numeric_limits<double>::infinity(),
       std::numeric_limits<double>::infinity(),S,I,C,N);
    REQUIRE((S.array() > 1e-6*diag).count() == 0);
    // Does not transversally (non-coplanar) pierce the input.
    Eigen::MatrixXi IF; Eigen::Array<bool,Eigen::Dynamic,1> CP;
    igl::find_intersections(CV,CF,V,F,/*first_only=*/false,IF,CP);
    int transversal = 0;
    for(int i = 0;i<CP.size();i++){ if(!CP(i)){ transversal++; } }
    REQUIRE(transversal == 0);
  }
}

TEST_CASE("lazy_cage: encloses input", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(test_common::data_path("cube.obj"),V,F);
  const double diag = igl::bounding_box_diagonal(V);

  const int num_faces = 40;
  const int grid_size = 24;
  Eigen::MatrixXd CV;
  Eigen::MatrixXi CF;
  double sigma = -1;
  // A cube must round out enough to be decimated this aggressively, so allow a
  // generous offset range.
  const bool ok = igl::lazy_cage(
    V,F,num_faces,grid_size,0.5*diag,8,/*use_qslim=*/false,
    igl::LAZY_CAGE_METRIC_SIGMA,igl::LAZY_CAGE_GRID_DENSE,igl::LAZY_CAGE_DISTANCE_SIGNED,CV,CF,sigma);

  REQUIRE(ok);
  REQUIRE(CF.rows() > 0);
  REQUIRE(CF.rows() <= num_faces);
  REQUIRE(sigma > 0.0);
  require_valid_cage(V,F,CV,CF);
}

// Regression: a coarse isosurfacing grid combined with a large target face
// count used to let the offset σ shrink below the grid resolution, where the
// dense offset already pierces the input; lazy_cage would then return a cage
// that intersected the input. The bisection must only accept σ that yield a
// genuinely valid (enclosing, non-piercing) cage.
TEST_CASE("lazy_cage: coarse grid does not pierce input", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(test_common::data_path("bunny_small.off"),V,F);

  const int num_faces = 800;
  const int grid_size = 48;   // deliberately coarse relative to the target
  Eigen::MatrixXd CV;
  Eigen::MatrixXi CF;
  double sigma = -1;
  const bool ok = igl::lazy_cage(V,F,num_faces,grid_size,CV,CF,sigma);

  REQUIRE(ok);
  REQUIRE(CF.rows() == num_faces);
  require_valid_cage(V,F,CV,CF);
}

TEST_CASE("lazy_cage: qslim variant is valid", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(test_common::data_path("cube.obj"),V,F);
  const double diag = igl::bounding_box_diagonal(V);

  Eigen::MatrixXd CV;
  Eigen::MatrixXi CF;
  double sigma = -1;
  const bool ok = igl::lazy_cage(
    V,F,/*num_faces=*/40,/*grid=*/24,/*max_sigma=*/0.5*diag,/*iters=*/8,
    /*use_qslim=*/true,igl::LAZY_CAGE_METRIC_SIGMA,igl::LAZY_CAGE_GRID_DENSE,igl::LAZY_CAGE_DISTANCE_SIGNED,CV,CF,sigma);
  REQUIRE(ok);
  require_valid_cage(V,F,CV,CF);
}

// The volume metric may select a larger σ than the σ metric if that yields a
// smaller-volume cage; either way the result must be a valid cage no larger in
// volume than the minimum-σ cage.
TEST_CASE("lazy_cage: volume metric is valid and no worse", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(test_common::data_path("bunny_small.off"),V,F);
  const double diag = igl::bounding_box_diagonal(V);
  const int num_faces = 800;
  const int grid_size = 48;

  Eigen::MatrixXd CVs, CVv;
  Eigen::MatrixXi CFs, CFv;
  double sig_s=-1, sig_v=-1;
  REQUIRE(igl::lazy_cage(
    V,F,num_faces,grid_size,0.1*diag,10,false,
    igl::LAZY_CAGE_METRIC_SIGMA,igl::LAZY_CAGE_GRID_DENSE,igl::LAZY_CAGE_DISTANCE_SIGNED,CVs,CFs,sig_s));
  REQUIRE(igl::lazy_cage(
    V,F,num_faces,grid_size,0.1*diag,10,false,
    igl::LAZY_CAGE_METRIC_VOLUME,igl::LAZY_CAGE_GRID_DENSE,igl::LAZY_CAGE_DISTANCE_SIGNED,CVv,CFv,sig_v));

  require_valid_cage(V,F,CVv,CFv);

  const auto volume = [](const Eigen::MatrixXd & CV, const Eigen::MatrixXi & CF)
  {
    Eigen::RowVector3d c; double vol=0; igl::centroid(CV,CF,c,vol);
    return std::abs(vol);
  };
  // The volume-optimal cage is at least as good (up to tiny tolerance) as the
  // minimum-σ cage under the volume objective.
  REQUIRE(volume(CVv,CFv) <= volume(CVs,CFs)*(1.0+1e-9));
}

// The sparse (Lipschitz octree) grid mode must produce a valid cage just like
// the dense mode, with a comparable offset.
TEST_CASE("lazy_cage: sparse grid mode is valid", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(test_common::data_path("cube.obj"),V,F);
  const double diag = igl::bounding_box_diagonal(V);

  Eigen::MatrixXd CVd, CVs;
  Eigen::MatrixXi CFd, CFs;
  double sig_d=-1, sig_s=-1;
  REQUIRE(igl::lazy_cage(
    V,F,40,32,0.5*diag,8,false,
    igl::LAZY_CAGE_METRIC_SIGMA,igl::LAZY_CAGE_GRID_DENSE,igl::LAZY_CAGE_DISTANCE_SIGNED,CVd,CFd,sig_d));
  REQUIRE(igl::lazy_cage(
    V,F,40,32,0.5*diag,8,false,
    igl::LAZY_CAGE_METRIC_SIGMA,igl::LAZY_CAGE_GRID_SPARSE,igl::LAZY_CAGE_DISTANCE_SIGNED,CVs,CFs,sig_s));

  require_valid_cage(V,F,CVd,CFd);
  require_valid_cage(V,F,CVs,CFs);
  // Same effective resolution, so the found offsets should be close (they can
  // differ because the two backends sample the field at different points).
  REQUIRE(sig_s == Approx(sig_d).epsilon(0.5));
}

// grid_size <= 0 selects the heuristic resolution and still yields a valid cage.
TEST_CASE("lazy_cage: automatic grid size", "[igl]")
{
  REQUIRE(igl::lazy_cage_default_grid_size(100) >= 24);
  REQUIRE(igl::lazy_cage_default_grid_size(100) <=
          igl::lazy_cage_default_grid_size(1600));
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(test_common::data_path("bunny_small.off"),V,F);
  Eigen::MatrixXd CV;
  Eigen::MatrixXi CF;
  double sigma = -1;
  REQUIRE(igl::lazy_cage(V,F,/*num_faces=*/200,/*grid_size=*/0,CV,CF,sigma));
  REQUIRE(CF.rows() == 200);
  require_valid_cage(V,F,CV,CF);
}

// Unsigned distance wraps an open / non-orientable input (here a flat sheet — a
// cloth-like mid-surface) in a closed cage, where signed distance is undefined.
TEST_CASE("lazy_cage: unsigned distance wraps an open sheet", "[igl]")
{
  // A flat n×n triangulated sheet on z=0 with an open boundary.
  const int n = 7;
  Eigen::MatrixXd V(n*n,3);
  for(int i=0;i<n;i++) for(int j=0;j<n;j++)
  {
    V.row(i*n+j) << double(i)/(n-1), double(j)/(n-1), 0.0;
  }
  Eigen::MatrixXi F(2*(n-1)*(n-1),3);
  int f=0;
  for(int i=0;i<n-1;i++) for(int j=0;j<n-1;j++)
  {
    const int a=i*n+j, b=i*n+j+1, c=(i+1)*n+j, d=(i+1)*n+j+1;
    F.row(f++)<<a,b,d;
    F.row(f++)<<a,d,c;
  }
  const double diag = igl::bounding_box_diagonal(V);

  Eigen::MatrixXd CV;
  Eigen::MatrixXi CF;
  double sigma = -1;
  const bool ok = igl::lazy_cage(
    V,F,/*num_faces=*/40,/*grid=*/32,/*max_sigma=*/0.4*diag,/*iters=*/9,
    /*use_qslim=*/false,igl::LAZY_CAGE_METRIC_SIGMA,igl::LAZY_CAGE_GRID_DENSE,
    igl::LAZY_CAGE_DISTANCE_UNSIGNED,CV,CF,sigma);
  REQUIRE(ok);
  REQUIRE(sigma > 0.0);
  require_valid_cage(V,F,CV,CF);

  // Sparse + unsigned should also work (and needs no winding number at all).
  Eigen::MatrixXd CVs;
  Eigen::MatrixXi CFs;
  double sig_s = -1;
  REQUIRE(igl::lazy_cage(
    V,F,40,32,0.4*diag,9,false,igl::LAZY_CAGE_METRIC_SIGMA,
    igl::LAZY_CAGE_GRID_SPARSE,igl::LAZY_CAGE_DISTANCE_UNSIGNED,CVs,CFs,sig_s));
  require_valid_cage(V,F,CVs,CFs);
}

#include <test_common.h>
#include <igl/lazy_cage.h>
#include <igl/read_triangle_mesh.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/signed_distance.h>
#include <igl/AABB.h>
#include <igl/tri_tri_intersect.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <limits>
#include <vector>

namespace
{
  // Count triangle pairs where a face of (AV,AF) transversally (non-coplanar)
  // intersects (BV,BF). Ground-truth-ish check using only core tools.
  int transversal_intersections(
    const Eigen::MatrixXd & AV, const Eigen::MatrixXi & AF,
    const Eigen::MatrixXd & BV, const Eigen::MatrixXi & BF)
  {
    igl::AABB<Eigen::MatrixXd,3> tree;
    tree.init(BV,BF);
    int count = 0;
    for(int f = 0; f < AF.rows(); f++)
    {
      const Eigen::RowVector3d A0 = AV.row(AF(f,0));
      const Eigen::RowVector3d A1 = AV.row(AF(f,1));
      const Eigen::RowVector3d A2 = AV.row(AF(f,2));
      Eigen::AlignedBox<double,3> box;
      box.extend(A0.transpose()).extend(A1.transpose()).extend(A2.transpose());
      std::vector<const igl::AABB<Eigen::MatrixXd,3>*> cand;
      tree.append_intersecting_leaves(box,cand);
      for(const auto * c : cand)
      {
        const int g = c->m_primitive;
        const Eigen::RowVector3d B0 = BV.row(BF(g,0));
        const Eigen::RowVector3d B1 = BV.row(BF(g,1));
        const Eigen::RowVector3d B2 = BV.row(BF(g,2));
        bool coplanar = false;
        Eigen::RowVector3d s,t;
        if(igl::tri_tri_intersection_test_3d(A0,A1,A2,B0,B1,B2,coplanar,s,t) &&
           !coplanar)
        {
          count++;
        }
      }
    }
    return count;
  }

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
    // Does not pierce the input.
    REQUIRE(transversal_intersections(CV,CF,V,F) == 0);
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
    V,F,num_faces,grid_size,0.5*diag,8,/*use_qslim=*/false,CV,CF,sigma);

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
    /*use_qslim=*/true,CV,CF,sigma);
  REQUIRE(ok);
  require_valid_cage(V,F,CV,CF);
}

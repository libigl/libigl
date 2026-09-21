#include <test_common.h>
#include <igl/triangle/remesh_at_points.h>
#include <igl/unique_edge_map.h>
#include <igl/boundary_facets.h>
#include <igl/is_edge_manifold.h>
#include <igl/placeholders.h>
#include <set>

namespace
{
  // Run igl::triangle::remesh_at_points on a 3D mesh and check all the
  // invariants promised by the header: VV starts with V, indices are in range,
  // the sample points appear as vertices at the right positions and are used by
  // the output mesh, and each original face is exactly tiled by its (positively
  // oriented) descendants.
  void check_remesh_at_points(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & B,
    const Eigen::VectorXi & FI,
    Eigen::MatrixXd & VV,
    Eigen::MatrixXi & FF,
    Eigen::VectorXi & J,
    Eigen::VectorXi & K)
  {
    REQUIRE(V.cols() == 3);
    igl::triangle::remesh_at_points(V,F,B,FI,VV,FF,J,K);

    // V is a prefix of VV
    REQUIRE(VV.cols() == V.cols());
    REQUIRE(VV.rows() >= V.rows());
    REQUIRE(VV.rows() <= V.rows() + B.rows());
    test_common::assert_eq(VV.topRows(V.rows()),V);

    // Sizes and index ranges
    REQUIRE(FF.cols() == 3);
    REQUIRE(FF.rows() > 0);
    REQUIRE(J.size() == FF.rows());
    REQUIRE(K.size() == B.rows());
    REQUIRE(FF.minCoeff() >= 0);
    REQUIRE(FF.maxCoeff() < VV.rows());
    REQUIRE(J.minCoeff() >= 0);
    REQUIRE(J.maxCoeff() < F.rows());
    REQUIRE(K.minCoeff() >= 0);
    REQUIRE(K.maxCoeff() < VV.rows());

    std::set<int> used(FF.data(),FF.data()+FF.size());
    for(int i = 0;i<B.rows();i++)
    {
      INFO("sample "<<i);
      // K(i) sits exactly where the barycentric coordinates say it should
      const Eigen::RowVector3d p =
        B(i,0)*V.row(F(FI(i),0)) +
        B(i,1)*V.row(F(FI(i),1)) +
        B(i,2)*V.row(F(FI(i),2));
      REQUIRE((VV.row(K(i))-p).norm() == Approx(0).margin(1e-14));
      // and it is actually referenced by the output mesh
      REQUIRE(used.count(K(i)) == 1);
    }

    // Every original face is exactly tiled by the faces that map to it, each of
    // which agrees with the original orientation.
    Eigen::VectorXd A(F.rows());
    for(int f = 0;f<F.rows();f++)
    {
      const Eigen::RowVector3d e1 = V.row(F(f,1))-V.row(F(f,0));
      const Eigen::RowVector3d e2 = V.row(F(f,2))-V.row(F(f,0));
      A(f) = e1.cross(e2).norm();
      REQUIRE(A(f) > 0);
    }
    Eigen::VectorXd A_acc = Eigen::VectorXd::Zero(F.rows());
    for(int i = 0;i<FF.rows();i++)
    {
      INFO("output face "<<i);
      const int f = J(i);
      Eigen::RowVector3d n =
        Eigen::RowVector3d(V.row(F(f,1))-V.row(F(f,0))).cross(
        Eigen::RowVector3d(V.row(F(f,2))-V.row(F(f,0))));
      n /= n.norm();
      const Eigen::RowVector3d m =
        Eigen::RowVector3d(VV.row(FF(i,1))-VV.row(FF(i,0))).cross(
        Eigen::RowVector3d(VV.row(FF(i,2))-VV.row(FF(i,0))));
      // Coplanar with the original face and consistently oriented (this also
      // rules out degenerate and overlapping output triangles).
      REQUIRE(m.norm() > 0);
      REQUIRE(m.dot(n) == Approx(m.norm()).margin(1e-12*A(f)));
      A_acc(f) += m.dot(n);
    }
    for(int f = 0;f<F.rows();f++)
    {
      INFO("original face "<<f);
      REQUIRE(A_acc(f) == Approx(A(f)).margin(1e-12*A.maxCoeff()));
    }
  }

  // A single triangle in 3D (not axis aligned, so the 2D barycentric detour is
  // actually exercised).
  void single_triangle(Eigen::MatrixXd & V, Eigen::MatrixXi & F)
  {
    V.resize(3,3);
    V <<
      0,0,0,
      2,0,1,
      0,3,2;
    F.resize(1,3);
    F<<0,1,2;
  }
}

TEST_CASE("triangle_remesh_at_points: face point", "[igl][triangle]")
{
  Eigen::MatrixXd V; Eigen::MatrixXi F;
  single_triangle(V,F);
  Eigen::MatrixXd B(1,3);
  B<<1.0/3.0,1.0/3.0,1.0/3.0;
  Eigen::VectorXi FI(1);
  FI<<0;

  Eigen::MatrixXd VV; Eigen::MatrixXi FF; Eigen::VectorXi J,K;
  check_remesh_at_points(V,F,B,FI,VV,FF,J,K);

  // Barycenter splits the triangle into three
  REQUIRE(VV.rows() == 4);
  REQUIRE(FF.rows() == 3);
  REQUIRE(K(0) == 3);
  REQUIRE((J.array() == 0).all());
}

TEST_CASE("triangle_remesh_at_points: vertex point", "[igl][triangle]")
{
  Eigen::MatrixXd V; Eigen::MatrixXi F;
  single_triangle(V,F);
  Eigen::MatrixXd B(1,3);
  B<<0,1,0;
  Eigen::VectorXi FI(1);
  FI<<0;

  Eigen::MatrixXd VV; Eigen::MatrixXi FF; Eigen::VectorXi J,K;
  check_remesh_at_points(V,F,B,FI,VV,FF,J,K);

  // Vertex-points are not inserted, just tracked in K
  REQUIRE(VV.rows() == V.rows());
  REQUIRE(FF.rows() == 1);
  REQUIRE(K(0) == F(0,1));
  test_common::assert_eq(FF,F);
  REQUIRE(J(0) == 0);
}

TEST_CASE("triangle_remesh_at_points: edge point on boundary edge", "[igl][triangle]")
{
  Eigen::MatrixXd V; Eigen::MatrixXi F;
  single_triangle(V,F);
  Eigen::MatrixXd B(1,3);
  B<<0,0.25,0.75;
  Eigen::VectorXi FI(1);
  FI<<0;

  Eigen::MatrixXd VV; Eigen::MatrixXi FF; Eigen::VectorXi J,K;
  check_remesh_at_points(V,F,B,FI,VV,FF,J,K);

  REQUIRE(VV.rows() == 4);
  REQUIRE(FF.rows() == 2);
  REQUIRE(K(0) == 3);
}

TEST_CASE("triangle_remesh_at_points: edge point splits both incident faces", "[igl][triangle]")
{
  // Two triangles sharing the diagonal 1-2
  Eigen::MatrixXd V(4,3);
  V <<
    0,0,0,
    1,0,0,
    0,1,0,
    1,1,0;
  Eigen::MatrixXi F(2,3);
  F<<
    0,1,2,
    1,3,2;
  Eigen::MatrixXd B(1,3);
  // Point on the shared edge (1,2) of face 0, 1/4 of the way from 1 to 2
  B<<0,0.75,0.25;
  Eigen::VectorXi FI(1);
  FI<<0;

  Eigen::MatrixXd VV; Eigen::MatrixXi FF; Eigen::VectorXi J,K;
  check_remesh_at_points(V,F,B,FI,VV,FF,J,K);

  REQUIRE(VV.rows() == 5);
  REQUIRE(FF.rows() == 4);
  // Both faces got split
  REQUIRE((J.array()==0).count() == 2);
  REQUIRE((J.array()==1).count() == 2);
  // Output is still manifold and, since the split edge was interior, still has
  // the same boundary
  REQUIRE(igl::is_edge_manifold(FF));
  Eigen::MatrixXi Ein,Eout;
  igl::boundary_facets(F,Ein);
  igl::boundary_facets(FF,Eout);
  REQUIRE(Eout.rows() == Ein.rows());
}

TEST_CASE("triangle_remesh_at_points: non-manifold edge", "[igl][triangle]")
{
  // Three triangles sharing the edge (0,1); note face 1 traverses the shared
  // edge in the opposite direction from faces 0 and 2.
  Eigen::MatrixXd V(5,3);
  V <<
    0, 0,0,
    1, 0,0,
    0, 1,0,
    0,-1,0,
    0, 0,1;
  Eigen::MatrixXi F(3,3);
  F<<
    0,1,2,
    1,0,3,
    0,1,4;

  SECTION("one point")
  {
    Eigen::MatrixXd B(1,3);
    B<<0.25,0.75,0;
    Eigen::VectorXi FI(1);
    FI<<0;

    Eigen::MatrixXd VV; Eigen::MatrixXi FF; Eigen::VectorXi J,K;
    check_remesh_at_points(V,F,B,FI,VV,FF,J,K);

    // The edge is split on all three incident faces
    REQUIRE(VV.rows() == 6);
    REQUIRE(FF.rows() == 6);
    for(int f = 0;f<3;f++) { REQUIRE((J.array()==f).count() == 2); }
  }

  SECTION("two points")
  {
    // Two points on the same shared edge: their order along the edge must be
    // consistent in every incident face, including the flipped one.
    Eigen::MatrixXd B(2,3);
    B<<
      0.75,0.25,0,
      0.25,0.75,0;
    Eigen::VectorXi FI(2);
    FI<<0,0;

    Eigen::MatrixXd VV; Eigen::MatrixXi FF; Eigen::VectorXi J,K;
    check_remesh_at_points(V,F,B,FI,VV,FF,J,K);

    REQUIRE(VV.rows() == 7);
    REQUIRE(FF.rows() == 9);
    for(int f = 0;f<3;f++) { REQUIRE((J.array()==f).count() == 3); }
  }
}

TEST_CASE("triangle_remesh_at_points: all edge midpoints of a cube", "[igl][triangle]")
{
  Eigen::MatrixXd V; Eigen::MatrixXi F;
  igl::read_triangle_mesh(test_common::data_path("cube.obj"),V,F);
  REQUIRE(igl::is_edge_manifold(F));

  Eigen::MatrixXi E,uE; Eigen::VectorXi EMAP,uEC,uEE;
  igl::unique_edge_map(F,E,uE,EMAP,uEC,uEE);

  // One point per unique edge, at its midpoint, expressed in the first
  // incident face we happen to see.
  std::vector<bool> seen(uE.rows(),false);
  std::vector<Eigen::RowVector3d> vB;
  std::vector<int> vFI;
  for(int c = 0;c<3;c++)
  {
    for(int f = 0;f<F.rows();f++)
    {
      const int ue = EMAP(c*F.rows()+f);
      if(seen[ue]) { continue; }
      seen[ue] = true;
      Eigen::RowVector3d b(0,0,0);
      b((c+1)%3) = 0.5;
      b((c+2)%3) = 0.5;
      vB.push_back(b);
      vFI.push_back(f);
    }
  }
  REQUIRE((int)vB.size() == uE.rows());
  Eigen::MatrixXd B(vB.size(),3);
  Eigen::VectorXi FI(vFI.size());
  for(int i = 0;i<(int)vB.size();i++) { B.row(i) = vB[i]; FI(i) = vFI[i]; }

  Eigen::MatrixXd VV; Eigen::MatrixXi FF; Eigen::VectorXi J,K;
  check_remesh_at_points(V,F,B,FI,VV,FF,J,K);

  // This is exactly one step of 4-1 subdivision
  REQUIRE(VV.rows() == V.rows()+uE.rows());
  REQUIRE(FF.rows() == 4*F.rows());
  for(int f = 0;f<F.rows();f++) { REQUIRE((J.array()==f).count() == 4); }
  // Still a closed manifold
  REQUIRE(igl::is_edge_manifold(FF));
  Eigen::MatrixXi Eout;
  igl::boundary_facets(FF,Eout);
  REQUIRE(Eout.rows() == 0);
}

TEST_CASE("triangle_remesh_at_points: many points per face", "[igl][triangle]")
{
  Eigen::MatrixXd V; Eigen::MatrixXi F;
  igl::read_triangle_mesh(test_common::data_path("decimated-knight.obj"),V,F);

  // A fixed table of strictly interior barycentric coordinates. All entries are
  // dyadic rationals so that each row sums to exactly 1 in floating point.
  const Eigen::Matrix<double,7,3> table = (Eigen::Matrix<double,7,3>() <<
    8.0/16.0, 4.0/16.0, 4.0/16.0,
    4.0/16.0, 8.0/16.0, 4.0/16.0,
    4.0/16.0, 4.0/16.0, 8.0/16.0,
   10.0/16.0, 4.0/16.0, 2.0/16.0,
    2.0/16.0,10.0/16.0, 4.0/16.0,
    4.0/16.0, 2.0/16.0,10.0/16.0,
    1.0/16.0, 7.0/16.0, 8.0/16.0).finished();
  REQUIRE((table.rowwise().sum().array() == 1.0).all());

  // Deterministically hand face f between 1 and 4 distinct points from the
  // table, at a face-dependent offset so faces get different configurations.
  std::vector<Eigen::RowVector3d> vB;
  std::vector<int> vFI;
  int max_per_face = 0;
  for(int f = 0;f<F.rows();f++)
  {
    const int count = (f % 4) + 1;
    max_per_face = std::max(max_per_face,count);
    for(int j = 0;j<count;j++)
    {
      vB.push_back(table.row((3*f + j) % table.rows()));
      vFI.push_back(f);
    }
  }
  REQUIRE(max_per_face == 4);
  Eigen::MatrixXd B(vB.size(),3);
  Eigen::VectorXi FI(vFI.size());
  for(int i = 0;i<(int)vB.size();i++) { B.row(i) = vB[i]; FI(i) = vFI[i]; }

  Eigen::MatrixXd VV; Eigen::MatrixXi FF; Eigen::VectorXi J,K;
  check_remesh_at_points(V,F,B,FI,VV,FF,J,K);

  // Interior points, so every one of them is a new vertex
  REQUIRE(VV.rows() == V.rows()+B.rows());
  // f interior points in a triangle means 2f+1 output triangles
  for(int f = 0;f<F.rows();f++)
  {
    INFO("original face "<<f);
    REQUIRE((J.array()==f).count() == 2*((f % 4) + 1) + 1);
  }
  REQUIRE(igl::is_edge_manifold(FF));
  Eigen::MatrixXi Eout;
  igl::boundary_facets(FF,Eout);
  REQUIRE(Eout.rows() == 0);
}

#include <test_common.h>
#include <igl/collapse_edge_would_intersect_mesh.h>
#include <igl/block_intersections_with_input.h>
#include <igl/block_self_intersections.h>
#include <igl/edge_flaps.h>
#include <igl/decimate.h>
#include <igl/qslim.h>
#include <igl/AABB.h>
#include <igl/decimate_callback_types.h>
#include <igl/read_triangle_mesh.h>
#include <algorithm>
#include <vector>

namespace
{
  // Unit octahedron centered at the origin.
  void octahedron(Eigen::MatrixXd & V, Eigen::MatrixXi & F)
  {
    V.resize(6,3);
    V <<
       1, 0, 0,
      -1, 0, 0,
       0, 1, 0,
       0,-1, 0,
       0, 0, 1,
       0, 0,-1;
    F.resize(8,3);
    F <<
      0,2,4,
      2,1,4,
      1,3,4,
      3,0,4,
      2,0,5,
      1,2,5,
      3,1,5,
      0,3,5;
  }

  int find_edge(const Eigen::MatrixXi & E, const int a, const int b)
  {
    int e = -1;
    const bool found =
      ( (E.array().col(0)==a && E.array().col(1)==b) ||
        (E.array().col(0)==b && E.array().col(1)==a) ).maxCoeff(&e);
    REQUIRE(found);
    return e;
  }
}

TEST_CASE(
  "collapse_edge_would_intersect_mesh: transversal vs distant",
  "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  octahedron(V,F);
  Eigen::VectorXi EMAP;
  Eigen::MatrixXi E,EF,EI;
  igl::edge_flaps(F,E,EMAP,EF,EI);

  const int e = find_edge(E,0,2);
  // Collapsing edge (0,2) places the merged vertex at the midpoint.
  const Eigen::RowVectorXd p = (0.5*(V.row(0)+V.row(2))).eval();

  // An obstacle triangle lying in the plane y = 0.167 that the reshaped
  // one-ring sweeps through (the fixed corners have y in {0,-1,1}, the merged
  // vertex has y=0.5, so faces incident on the moved corners cross this plane).
  Eigen::MatrixXd VB(3,3);
  VB <<
    -0.5, 0.167, -0.2,
     0.3, 0.167, -0.2,
    -0.1, 0.167,  0.9;
  Eigen::MatrixXi FB(1,3);
  FB << 0,1,2;
  igl::AABB<Eigen::MatrixXd,3> tree;
  tree.init(VB,FB);
  REQUIRE(igl::collapse_edge_would_intersect_mesh(
    e,p,V,F,E,EMAP,EF,EI,VB,FB,tree));

  // A distant obstacle (far above) is never intersected by the collapse.
  Eigen::MatrixXd VB2 = VB;
  VB2.col(1).array() += 5.0;
  igl::AABB<Eigen::MatrixXd,3> tree2;
  tree2.init(VB2,FB);
  REQUIRE_FALSE(igl::collapse_edge_would_intersect_mesh(
    e,p,V,F,E,EMAP,EF,EI,VB2,FB,tree2));
}

TEST_CASE(
  "collapse_edge_would_intersect_mesh: adjacent obstacle ignored",
  "[igl]")
{
  // An obstacle triangle that shares a fixed (un-moved) corner of a reshaped
  // one-ring triangle only touches it tangentially and must not be reported.
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  octahedron(V,F);
  Eigen::VectorXi EMAP;
  Eigen::MatrixXi E,EF,EI;
  igl::edge_flaps(F,E,EMAP,EF,EI);
  const int e = find_edge(E,0,2);
  const Eigen::RowVectorXd p = (0.5*(V.row(0)+V.row(2))).eval();
  // Obstacle triangle incident on octahedron vertex 4=(0,0,1), which is a fixed
  // corner of several reshaped one-ring faces. It fans away from the mesh so it
  // only touches at that shared vertex.
  Eigen::MatrixXd VB(3,3);
  VB <<
    0,0,1,
    3,1,4,
    1,3,4;
  Eigen::MatrixXi FB(1,3);
  FB << 0,1,2;
  igl::AABB<Eigen::MatrixXd,3> tree;
  tree.init(VB,FB);
  REQUIRE_FALSE(igl::collapse_edge_would_intersect_mesh(
    e,p,V,F,E,EMAP,EF,EI,VB,FB,tree));
}

TEST_CASE("block decorators: qslim end-to-end", "[igl]")
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  // A closed manifold mesh.
  igl::read_triangle_mesh(test_common::data_path("cube.obj"),V,F);

  const int target = std::max<int>(8,F.rows()/2);

  // Build an outward, disjoint "shell" obstacle by scaling the cube about its
  // centroid. The decimating cube lives strictly inside it and should never be
  // blocked by it, so results should match plain qslim.
  const Eigen::RowVector3d c = V.colwise().mean();
  Eigen::MatrixXd VB = (V.rowwise()-c)*1.5;
  VB.rowwise() += c;
  Eigen::MatrixXi FB = F;

  std::vector<igl::decimate_pre_post_collapse_callbacks_decorator> decorators;
  decorators.push_back(igl::block_self_intersections());
  decorators.push_back(igl::block_intersections_with_input(VB,FB));

  Eigen::MatrixXd U;
  Eigen::MatrixXi G;
  Eigen::VectorXi J,I;
  igl::qslim(V,F,target,decorators,U,G,J,I);
  REQUIRE(G.rows() > 0);
  REQUIRE(G.rows() <= F.rows());
  // Indices are in bounds.
  REQUIRE(G.maxCoeff() < U.rows());
  REQUIRE(G.minCoeff() >= 0);

  // The empty-decorator overload should behave like plain qslim.
  std::vector<igl::decimate_pre_post_collapse_callbacks_decorator> none;
  Eigen::MatrixXd U2; Eigen::MatrixXi G2; Eigen::VectorXi J2,I2;
  igl::qslim(V,F,target,none,U2,G2,J2,I2);
  REQUIRE(G2.rows() > 0);
}

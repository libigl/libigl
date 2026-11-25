#include <test_common.h>
#include <igl/eytzinger_aabb.h>
#include <igl/upsample.h>
#include <cmath>
#include <limits>
#include <vector>

namespace
{
  // Number of primitives stored below node i. Also checks that every internal
  // node has its children inside the arrays, splits its primitives evenly, and
  // stores the union of its children's boxes.
  int check_subtree(
    const Eigen::MatrixXd & B1,
    const Eigen::MatrixXd & B2,
    const Eigen::VectorXi & leaf,
    const int i)
  {
    if(leaf(i) == -2) { return 0; }
    if(leaf(i) >= 0) { return 1; }
    REQUIRE(leaf(i) == -1);
    const int left = 2*i + 1;
    const int right = 2*i + 2;
    REQUIRE(right < leaf.size());
    const int num_left = check_subtree(B1,B2,leaf,left);
    const int num_right = check_subtree(B1,B2,leaf,right);
    REQUIRE(num_left > 0);
    REQUIRE(num_right > 0);
    REQUIRE(std::abs(num_left - num_right) <= 1);
    test_common::assert_eq(B1.row(i),B1.row(left).cwiseMin(B1.row(right)));
    test_common::assert_eq(B2.row(i),B2.row(left).cwiseMax(B2.row(right)));
    return num_left + num_right;
  }

  // Build a tree over the boxes (PB1,PB2) and check that it is complete.
  void check_tree(const Eigen::MatrixXd & PB1, const Eigen::MatrixXd & PB2)
  {
    Eigen::MatrixXd B1,B2;
    Eigen::VectorXi leaf;
    igl::eytzinger_aabb(PB1,PB2,B1,B2,leaf);
    const int m = PB1.rows();
    const int num_nodes = 1 << ((int)std::ceil(std::log2(m))+1);
    REQUIRE(leaf.size() == num_nodes);
    REQUIRE(B1.rows() == num_nodes);
    REQUIRE(B2.rows() == num_nodes);
    // Every primitive is in exactly one leaf, and that leaf stores its box.
    std::vector<int> count(m,0);
    for(int i = 0; i < leaf.size(); i++)
    {
      if(leaf(i) < 0) { continue; }
      REQUIRE(leaf(i) < m);
      count[leaf(i)]++;
      test_common::assert_eq(B1.row(i),PB1.row(leaf(i)));
      test_common::assert_eq(B2.row(i),PB2.row(leaf(i)));
    }
    for(int p = 0; p < m; p++) { REQUIRE(count[p] == 1); }
    REQUIRE(check_subtree(B1,B2,leaf,0) == m);
  }
}

TEST_CASE("eytzinger_aabb: tied_centers_sphere", "[igl]")
{
  // An octahedron subdivided onto the unit sphere is symmetric, so many face
  // boxes share a center coordinate.
  Eigen::MatrixXd V(6,3);
  V<<1,0,0, -1,0,0, 0,1,0, 0,-1,0, 0,0,1, 0,0,-1;
  Eigen::MatrixXi F(8,3);
  F<<0,2,4, 2,1,4, 1,3,4, 3,0,4, 2,0,5, 1,2,5, 3,1,5, 0,3,5;
  for(int level = 0; level < 5; level++)
  {
    Eigen::MatrixXd NV;
    Eigen::MatrixXi NF;
    igl::upsample(V,F,NV,NF);
    V = NV;
    F = NF;
    V.rowwise().normalize();
  }
  REQUIRE(F.rows() == 8192);
  Eigen::MatrixXd PB1(F.rows(),3),PB2(F.rows(),3);
  for(int f = 0; f < F.rows(); f++)
  {
    PB1.row(f) = V.row(F(f,0)).cwiseMin(V.row(F(f,1))).cwiseMin(V.row(F(f,2)));
    PB2.row(f) = V.row(F(f,0)).cwiseMax(V.row(F(f,1))).cwiseMax(V.row(F(f,2)));
  }
  check_tree(PB1,PB2);
}

TEST_CASE("eytzinger_aabb: identical_boxes", "[igl]")
{
  // Every center ties with every other, for sizes on both sides of a power
  // of two.
  for(const int m : {1,2,3,5,31,32,33})
  {
    const Eigen::MatrixXd PB1 = Eigen::MatrixXd::Constant(m,3,-1);
    const Eigen::MatrixXd PB2 = Eigen::MatrixXd::Constant(m,3,2);
    check_tree(PB1,PB2);
  }
}

TEST_CASE("eytzinger_aabb: grid_2d", "[igl]")
{
  // Unit cells of a 7 by 9 grid: each center coordinate repeats in a whole row
  // or column.
  const int nx = 7;
  const int ny = 9;
  Eigen::MatrixXd PB1(nx*ny,2),PB2(nx*ny,2);
  for(int x = 0; x < nx; x++)
  {
    for(int y = 0; y < ny; y++)
    {
      PB1.row(x*ny+y) << x,y;
      PB2.row(x*ny+y) << x+1,y+1;
    }
  }
  check_tree(PB1,PB2);
}

TEST_CASE("eytzinger_aabb: empty", "[igl]")
{
  Eigen::MatrixXd B1,B2;
  Eigen::VectorXi leaf;
  igl::eytzinger_aabb(Eigen::MatrixXd(0,3),Eigen::MatrixXd(0,3),B1,B2,leaf);
  REQUIRE(leaf.size() == 0);
  REQUIRE(B1.rows() == 0);
  REQUIRE(B2.rows() == 0);
}

TEST_CASE("eytzinger_aabb: nan_boxes_stay_in_bounds", "[igl]")
{
  // Boxes with NaN coordinates cannot be ordered meaningfully, but every
  // primitive must still land in exactly one leaf inside the arrays.
  const int m = 21;
  Eigen::MatrixXd PB1 = Eigen::MatrixXd::Zero(m,3);
  Eigen::MatrixXd PB2 = Eigen::MatrixXd::Ones(m,3);
  for(int p = 0; p < m; p++)
  {
    PB1(p,0) = p%4;
    PB2(p,0) = p%4 + 1;
  }
  for(const int p : {0,7,8,20})
  {
    PB1.row(p).setConstant(std::numeric_limits<double>::quiet_NaN());
    PB2.row(p).setConstant(std::numeric_limits<double>::quiet_NaN());
  }
  Eigen::MatrixXd B1,B2;
  Eigen::VectorXi leaf;
  IGL_PUSH_FPE;
  igl::eytzinger_aabb(PB1,PB2,B1,B2,leaf);
  IGL_POP_FPE;
  REQUIRE(leaf.size() == 64);
  std::vector<int> count(m,0);
  for(int i = 0; i < leaf.size(); i++)
  {
    if(leaf(i) >= 0)
    {
      REQUIRE(leaf(i) < m);
      count[leaf(i)]++;
    }
  }
  for(int p = 0; p < m; p++) { REQUIRE(count[p] == 1); }
}

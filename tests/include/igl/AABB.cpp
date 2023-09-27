#include <igl/AABB.h>
#include <test_common.h>

TEST_CASE("AABB: find_2d", "[igl]")
{
  Eigen::MatrixXd V(6,2);
  V << 
  0,0,
  1,0,
  0,1,
  2,1,
  2,2,
  1,2;
  Eigen::MatrixXi F(4,3);
  F<< 
  2,0,1,
  2,1,5,
  5,3,4,
  5,1,3;
  igl::AABB<Eigen::MatrixXd,2> tree;
  tree.init(V,F);
  Eigen::RowVector2d q(0.5,0.5);
  std::vector<int> r = tree.find(V,F,q);
  REQUIRE(r.size() == 2);
  REQUIRE(r[0] == 0);
  REQUIRE(r[1] == 1);
}

TEST_CASE("AABB: find_3d", "[igl]")
{
  Eigen::MatrixXd V(7,3);
  V << 0,0,1,
    1,0,1,
    0,1,1,
    2,1,1,
    2,2,1,
    1,2,1,
    0,0,0;

  Eigen::MatrixXi F(4,4);
  F << 
    0,1,2,6,
    1,3,2,6,
    3,4,5,6,
    3,5,1,6;

  igl::AABB<Eigen::MatrixXd,3> tree;
  tree.init(V,F);
  Eigen::RowVector3d q(0.5,0.5,1.0);
  std::vector<int> r = tree.find(V,F,q);
  REQUIRE(r.size() == 2);
  REQUIRE(r[0] == 0);
  REQUIRE(r[1] == 1);
}

#include <test_common.h>
#include <igl/find_intersections.h>
#include <Eigen/Core>

TEST_CASE("find_intersections: crossing vs disjoint", "[igl]")
{
  // A single triangle in the z=0 plane.
  Eigen::MatrixXd V1(3,3);
  V1 << 0,0,0,  2,0,0,  0,2,0;
  Eigen::MatrixXi F1(1,3); F1 << 0,1,2;

  // A triangle that passes transversally through it.
  Eigen::MatrixXd V2(3,3);
  V2 << 0.5,0.5,-1,  0.5,0.5,1,  1.5,0.2,0;
  Eigen::MatrixXi F2(1,3); F2 << 0,1,2;

  Eigen::MatrixXi IF; Eigen::Array<bool,Eigen::Dynamic,1> CP;
  REQUIRE(igl::find_intersections(V1,F1,V2,F2,/*first_only=*/false,IF,CP));
  REQUIRE(IF.rows() == 1);
  REQUIRE(CP(0) == false);
  REQUIRE(IF(0,0) == 0);
  REQUIRE(IF(0,1) == 0);

  // Move the second triangle far away: no intersection.
  Eigen::MatrixXd V2b = V2;
  V2b.col(2).array() += 10.0;
  Eigen::MatrixXi IFb; Eigen::Array<bool,Eigen::Dynamic,1> CPb;
  REQUIRE_FALSE(igl::find_intersections(V1,F1,V2b,F2,false,IFb,CPb));
  REQUIRE(IFb.rows() == 0);
}

TEST_CASE("find_intersections: first_only early-out", "[igl]")
{
  // Two coplanar unit squares (2 triangles each) offset in x so several pairs
  // of triangles overlap.
  Eigen::MatrixXd V1(4,3);
  V1 << 0,0,0, 1,0,0, 1,1,0, 0,1,0;
  Eigen::MatrixXi F1(2,3); F1 << 0,1,2,  0,2,3;
  Eigen::MatrixXd V2 = V1;
  V2.col(0).array() += 0.5;               // overlapping, coplanar
  Eigen::MatrixXi F2 = F1;

  Eigen::MatrixXi IF_all; Eigen::Array<bool,Eigen::Dynamic,1> CP_all;
  REQUIRE(igl::find_intersections(V1,F1,V2,F2,false,IF_all,CP_all));

  Eigen::MatrixXi IF_first; Eigen::Array<bool,Eigen::Dynamic,1> CP_first;
  REQUIRE(igl::find_intersections(V1,F1,V2,F2,true,IF_first,CP_first));
  // first_only returns exactly one pair; the full search returns at least that.
  REQUIRE(IF_first.rows() == 1);
  REQUIRE(IF_all.rows() >= IF_first.rows());
}

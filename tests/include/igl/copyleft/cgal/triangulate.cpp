#include <test_common.h>
#include <igl/copyleft/cgal/triangulate.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

TEST_CASE("igl_copyleft_cgal_triangulate: sqannulus", "[igl/copyleft/cgal]")
{
  Eigen::MatrixXd V(8,2);
  V<<0,0,3,0,3,3,0,3,
     1,1,1,2,2,2,2,1;
  Eigen::MatrixXi E(8,2);
  E<<0,1,1,2,2,3,3,0,
     4,5,5,6,6,7,7,4;
  Eigen::MatrixXd H(1,2);
  H<<1.5,1.5;
  Eigen::MatrixXd TV;
  Eigen::MatrixXi TF;
  igl::copyleft::cgal::triangulate<CGAL::Epeck>(V,E,H,false,TV,TF);
  Eigen::MatrixXd gt_TV(8,2);
  gt_TV<<
  0,0,
  3,0,
  3,3,
  0,3,
  1,1,
  1,2,
  2,2,
  2,1;
  Eigen::MatrixXi gt_TF(8,3);
  gt_TF<<
  7,4,0,
  7,0,1,
  3,0,4,
  3,5,6,
  3,4,5,
  2,6,1,
  2,3,6,
  6,7,1;
  test_common::assert_eq(TV,gt_TV);
  test_common::assert_eq(TF,gt_TF);
}


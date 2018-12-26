#include <test_common.h>
#include <igl/boundary_facets.h>
#include <igl/sort.h>
#include <igl/sortrows.h>
#include <igl/setxor.h>


TEST_CASE("boundary_facets: single_tet", "[igl]")
{
  const Eigen::MatrixXi T = 
    (Eigen::MatrixXi(1,4)<<0,1,2,3).finished();
  Eigen::MatrixXi F;
  Eigen::MatrixXi J;
  Eigen::MatrixXi K;
  igl::boundary_facets(T,F);
  REQUIRE( F.rows () == 4 );
  // sorted! (unoriented)
  const Eigen::MatrixXi sorted_Fgt = 
    (Eigen::MatrixXi(4,3) << 
     0,1,2,
     0,1,3,
     0,2,3,
     1,2,3).finished();
  Eigen::MatrixXi sorted_F;
  {
    Eigen::MatrixXi _1;
    igl::sort(F,2,true, sorted_F,_1);
    igl::sortrows(Eigen::MatrixXi(sorted_F),true,sorted_F,_1);
  }
  test_common::assert_eq(sorted_Fgt,sorted_F);
}

TEST_CASE("boundary_facets: single_cube", "[igl]")
{
  const Eigen::MatrixXi T = 
    (Eigen::MatrixXi(6,4)<<
    0,1,7,5,
    0,7,4,5,
    0,1,3,7,
    0,3,2,7,
    0,6,4,7,
    0,2,6,7).finished();
  const Eigen::MatrixXi Fgt = 
    (Eigen::MatrixXi(12,3)<<
    0,5,4,
    0,1,5,
    6,7,2,
    7,3,2,
    4,6,0,
    6,2,0,
    1,7,5,
    1,3,7,
    0,3,1,
    0,2,3,
    5,7,4,
    7,6,4).finished();
  Eigen::MatrixXi F;
  igl::boundary_facets(T,F);
  const auto sortF = [](const Eigen::MatrixXi & F)-> Eigen::MatrixXi
  {
    Eigen::MatrixXi sorted_F;
    Eigen::MatrixXi _1;
    igl::sort(F,2,true, sorted_F,_1);
    igl::sortrows(Eigen::MatrixXi(sorted_F),true,sorted_F,_1);
    return sorted_F;
  };
  Eigen::MatrixXi sorted_F = sortF(F);
  Eigen::MatrixXi sorted_Fgt = sortF(Fgt);
  test_common::assert_eq(sorted_Fgt,sorted_F);
}

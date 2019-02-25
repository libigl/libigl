#include <test_common.h>
#include <igl/boundary_facets.h>
#include <igl/sort.h>
#include <igl/sortrows.h>
#include <igl/setxor.h>

#include <igl/matlab_format.h>
#include <iostream>


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
  Eigen::VectorXi J,K;
  igl::boundary_facets(T,F,J,K);
  const auto sortF = [](const Eigen::MatrixXi & F)-> Eigen::MatrixXi
  {
    Eigen::MatrixXi sorted_F;
    igl::sort(F,2,true, sorted_F);
    igl::sortrows(Eigen::MatrixXi(sorted_F),true,sorted_F);
    return sorted_F;
  };
  Eigen::MatrixXi sorted_F = sortF(F);
  Eigen::MatrixXi sorted_Fgt = sortF(Fgt);
  test_common::assert_eq(sorted_Fgt,sorted_F);
  for(int f = 0;f<F.rows();f++)
  {
    Eigen::RowVector3i Ff;
    igl::sort(F.row(f),2,true,Ff);
    Eigen::RowVector3i Gf(
      T(J(f), (K(f)+1)%4 ),
      T(J(f), (K(f)+2)%4 ),
      T(J(f), (K(f)+3)%4 ));
    igl::sort(Eigen::RowVector3i(Gf),2,true,Gf);
    test_common::assert_eq(Ff,Gf);
  }
}

TEST_CASE("boundary_facets: non-manifold", "[igl]")
{
  const auto F = (Eigen::MatrixXi(3,3)<<0,1,2,1,0,3,0,1,4).finished();
  auto Egt = (Eigen::MatrixXi(6,2)<<1,2,2,0,3,1,0,3,4,0,1,4).finished();
  igl::sortrows(Eigen::MatrixXi(Egt),true,Egt);
  Eigen::MatrixXi E;
  igl::boundary_facets(F,E);
  std::cerr<<igl::matlab_format(E,"E")<<std::endl;
  igl::sortrows(Eigen::MatrixXi(E),true,E);
  test_common::assert_eq(Egt,E);
}

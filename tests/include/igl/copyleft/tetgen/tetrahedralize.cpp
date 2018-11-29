#include <test_common.h>

#include <igl/setxor.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>

TEST_CASE("tetrahedralize: single_tet", "[igl/copyleft/tetgen]")
{
  const Eigen::MatrixXd V = (Eigen::MatrixXd(4,3)<<
    0,0,0,
    1,0,0,
    0,1,0,
    0,0,1).finished();
  Eigen::MatrixXi F(0,3);
  Eigen::MatrixXd TV;
  Eigen::MatrixXi TT,TF;
  igl::copyleft::tetgen::tetrahedralize(V,F,"cpQ",TV,TT,TF);
  REQUIRE (4 == TV.rows());
  REQUIRE (1 == TT.rows());
  REQUIRE (4 == TF.rows());
  Eigen::MatrixXi TTcorrect = (Eigen::MatrixXi(1,4)<<0,1,2,3).finished();
  {
    Eigen::VectorXi TT_XOR,IA,IB;
    igl::setxor(TT,TTcorrect,TT_XOR,IA,IB);
    REQUIRE (0 == TT_XOR.size());
  }
  test_common::assert_eq(TV,V);
}

TEST_CASE("tetrahedralize: schoenhardt", "[igl/copyleft/tetgen]")
{
  const Eigen::MatrixXd V = (Eigen::MatrixXd(6,3)<<
    -246.2,-43.412,500,
     160.7,-191.51,500,
    85.505, 234.92,500,
       250,      0,  0,
      -125, 216.51,  0,
      -125,-216.51,  0).finished();
  const Eigen::MatrixXi F = (Eigen::MatrixXi(8,3)<<
    1,2,0,
    5,0,2,
    3,1,0,
    4,2,1,
    0,5,3,
    1,3,4,
    2,4,5,
    3,5,4).finished();
  Eigen::MatrixXd TV;
  Eigen::MatrixXi TT,TF;
  igl::copyleft::tetgen::tetrahedralize(V,F,"pQ",TV,TT,TF);
  REQUIRE (V.rows() <= TV.rows());
}

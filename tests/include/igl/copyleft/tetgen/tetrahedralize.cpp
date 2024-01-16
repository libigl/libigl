#include <test_common.h>

#include <igl/setxor.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/unique_simplices.h>
#include <igl/matlab_format.h>
#include <iostream>

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

TEST_CASE("tetrahedralize: two_tets", "[igl/copyleft/tetgen]")
{
  const Eigen::MatrixXd V = (Eigen::MatrixXd(5,3)<<
    2,1,0,
    1,-1,0,
    -1,-1,0,
    -1,1,0,
     0,0,1 ).finished();
  Eigen::MatrixXi F(0,3);
  Eigen::MatrixXd TV;
  Eigen::MatrixXi TT,TF;
  Eigen::MatrixXd H,R;
  Eigen::VectorXi FM,PT,TM,TR;
  Eigen::VectorXi VM = (Eigen::VectorXi(5)<< 0, 1, 2, 3, 4).finished();
  Eigen::MatrixXi TN,FT;
  int num_regions;
  igl::copyleft::tetgen::tetrahedralize(
    V,F,H,VM,FM,R,"cpnnfmQ",TV,TT,TF,TM,TR,TN,PT,FT,num_regions);

  REQUIRE (5 == TV.rows());
  REQUIRE (2 == TT.rows());
  REQUIRE (7 == TF.rows());
  REQUIRE (TV.rows() == TM.rows());
  REQUIRE (TN.rows() == TT.rows());
  REQUIRE (TN.cols() == 4);
  REQUIRE((TN.row(0).array() == 1).any());
  REQUIRE((TN.row(1).array() == 0).any());
  REQUIRE(FT.rows() == TF.rows());


  REQUIRE(PT.size() == TV.rows());
  REQUIRE(PT(0) == 0);
  REQUIRE(PT(2) == 1);

  Eigen::MatrixXi TTcorrect(2,4);
  TTcorrect<< 
    4,0,1,3,
    4,3,1,2;
  {
    Eigen::MatrixXi TT_combined = Eigen::MatrixXi(TT.rows()+TTcorrect.rows(),TT.cols());
    TT_combined<<TT,TTcorrect;
    Eigen::MatrixXi TT_unique;
    igl::unique_simplices(TT_combined,TT_unique);
    REQUIRE(TT_unique.rows() == TT.rows());
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

TEST_CASE("tetrahedralize: quad_face", "[igl/copyleft/tetgen]")
{
  const Eigen::MatrixXd V = (Eigen::MatrixXd(5,3)<<
    -1.0,0.0,0.0,
    0.0,-1.0,0.0,
    1.0,0.0,0.0,
    0.0,1.0,0.0,
    0.0,0.0,1.0).finished();
  const Eigen::MatrixXi F = (Eigen::MatrixXi(5,4)<<
    0,1,2,3,
    0,1,4,-1,
    1,2,4,-1,
    2,3,4,-1,
    3,0,4,-1).finished();
  Eigen::MatrixXd TV;
  Eigen::MatrixXi TT,TF;
  igl::copyleft::tetgen::tetrahedralize(V,F,"pQ",TV,TT,TF);
  REQUIRE (TT.rows() == 2);
}

TEST_CASE("tetrahedralize: embedded_segment", "[igl/copyleft/tetgen]")
{
  const Eigen::MatrixXd V = (Eigen::MatrixXd(6,3)<<
    -0.5,0.5,0.0,
    0.0,-1.0,0.0,
    0.5,0.5,0.0,
    0.0,0.0,1.0,
    0.0,0.0,0.2,
    0.0,0.0,0.9).finished();
  const Eigen::MatrixXi F = (Eigen::MatrixXi(5,3)<<
    0,1,2,
    1,3,2,
    0,2,3,
    1,0,3,
    4,5,-1).finished();
  Eigen::MatrixXd TV;
  Eigen::MatrixXi TT,TF;
  igl::copyleft::tetgen::tetrahedralize(V,F,"pQ",TV,TT,TF);
  REQUIRE (TT.rows() == 7);
}
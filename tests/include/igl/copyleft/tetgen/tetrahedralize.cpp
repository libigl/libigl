#include <test_common.h>

#include <igl/setxor.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>

TEST(tetrahedralize, single_tet) {
  const Eigen::MatrixXd V = (Eigen::MatrixXd(4,3)<<
    0,0,0,
    1,0,0,
    0,1,0,
    0,0,1).finished();
  Eigen::MatrixXi F(0,3);
  Eigen::MatrixXd TV;
  Eigen::MatrixXi TT,TF;
  igl::copyleft::tetgen::tetrahedralize(V,F,"cpQ",TV,TT,TF);
  ASSERT_EQ(TV.rows() , 4);
  ASSERT_EQ(TT.rows() , 1);
  ASSERT_EQ(TF.rows() , 4);
  Eigen::MatrixXi TTcorrect = (Eigen::MatrixXi(1,4)<<0,1,2,3).finished();
  {
    Eigen::VectorXi TT_XOR,IA,IB;
    igl::setxor(TT,TTcorrect,TT_XOR,IA,IB);
    ASSERT_EQ(TT_XOR.size(),0);
  }
  test_common::assert_eq(TV,V);
}

TEST(tetrahedralize, schoenhardt) {
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
  ASSERT_GE(TV.rows() , V.rows());
}

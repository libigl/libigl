#include <test_common.h>
#include <igl/split_nonmanifold.h>
#include <igl/matlab_format.h>
#include <iostream>

TEST_CASE("split_nonmanifold: fan", "[igl]")
{
  using namespace igl;
  Eigen::MatrixXd V(7,3);
  V << 0,0,0,
       1,0,0,
      -1,0,0,
       0,1,0,
       0,0,1,
       0,0,-1,
       1,0,1;
  Eigen::MatrixXi F(5,3);
  F<<0,1,3,
     0,3,2,
     0,4,3,
     0,3,5,
     0,3,6;
  Eigen::MatrixXd SV;
  Eigen::MatrixXi SF;
  Eigen::VectorXi SVI;
  igl::split_nonmanifold(V,F,SV,SF,SVI);
  Eigen::MatrixXd SVgt(13,3);
  SVgt<<
    0,0,0,
    0,0,0,
    0,0,0,
    0,0,0,
    1,0,0,
    0,1,0,
    0,0,1,
    0,1,0,
    0,1,0,
    -1,0,0,
    0,1,0,
    0,0,-1,
    1,0,1;
  Eigen::MatrixXi SFgt(5,3);
  SFgt<<
    0,4,5,
    0,5,9,
    1,6,10,
    2,7,11,
    3,8,12;
  Eigen::VectorXi SVIgt(13);
  SVIgt<<0,0,0,0,1,3,4,3,3,2,3,5,6;
  test_common::assert_eq(SV,SVgt);
  test_common::assert_eq(SF,SFgt);
  test_common::assert_eq(SVI,SVIgt);
}

TEST_CASE("split_nonmanifold: vertex", "[igl]")
{
  using namespace igl;
  Eigen::MatrixXd V(5,2);
  V << 0,0,
       1,0,
       0,1,
       2,0,
       1,1;
  Eigen::MatrixXi F(2,3);
  F<<0,1,2,
     1,3,4;
  Eigen::MatrixXd SV;
  Eigen::MatrixXi SF;
  Eigen::VectorXi SVI;
  igl::split_nonmanifold(V,F,SV,SF,SVI);
  Eigen::MatrixXd SVgt(6,2);
  SVgt<<
    0,0,
    1,0,
    1,0,
    2,0,
    0,1,
    1,1;
  Eigen::MatrixXi SFgt(2,3);
  SFgt<<
    0,2,4,
    1,3,5;
  Eigen::VectorXi SVIgt(6);
  SVIgt << 0, 1, 1, 3, 2, 4;
  test_common::assert_eq(SV,SVgt);
  test_common::assert_eq(SF,SFgt);
  test_common::assert_eq(SVI,SVIgt);
}


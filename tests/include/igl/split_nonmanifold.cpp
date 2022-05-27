#include <test_common.h>
#include <igl/split_nonmanifold.h>
#include <igl/matlab_format.h>
#include <iostream>

TEST_CASE("split_nonmanifold: edge-fan", "[igl]")
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

TEST_CASE("split_nonmanifold: vertex-boundary", "[igl]")
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

TEST_CASE("split_nonmanifold: edge-disk-flap", "[igl]")
{
  using namespace igl;
  Eigen::MatrixXd V(6,3);
  V<<
    0,0,0,
    1,0,0,
    0,1,0,
    -1,0,0,
    0,-1,0,
    0,0,1;
  Eigen::MatrixXi F(5,3);
  F<<
    0,1,2,
    0,2,3,
    0,3,4,
    0,4,1,
    0,5,1;
  Eigen::MatrixXd SV;
  Eigen::MatrixXi SF;
  Eigen::VectorXi SVI;
  igl::split_nonmanifold(V,F,SV,SF,SVI);
  Eigen::MatrixXd SVgt(8,3);
  SVgt<<
    0,0,0,
    0,0,0,
    1,0,0,
    0,1,0,
    -1,0,0,
    0,-1,0,
    0,0,1,
    1,0,0;
  Eigen::MatrixXi SFgt(5,3);
  SFgt<<
    0,2,3,
    0,3,4,
    0,4,5,
    0,5,2,
    1,6,7;
  Eigen::VectorXi SVIgt(8);
  SVIgt<<  0,  0,  1,  2,  3,  4,  5,  1;
  test_common::assert_eq(SV,SVgt);
  test_common::assert_eq(SF,SFgt);
  test_common::assert_eq(SVI,SVIgt);
}

TEST_CASE("split_nonmanifold: edge-disk-tent", "[igl]")
{
  using namespace igl;
  Eigen::MatrixXd V(5,3);
  V<<
    0,0,0,
    1,0,0,
    -1,1,0,
    0,-1,0,
    0,0,1;
  Eigen::MatrixXi F(5,3);
  F<<
    0,1,2,
    0,2,3,
    0,3,1,
    0,4,1,
    1,4,3;
  Eigen::MatrixXd SV;
  Eigen::MatrixXi SF;
  Eigen::VectorXi SVI;
  igl::split_nonmanifold(V,F,SV,SF,SVI);
  Eigen::MatrixXd SVgt(8,3);
  SVgt<<
    0,0,0,
    0,0,0,
    1,0,0,
    1,0,0,
    -1,1,0,
    0,-1,0,
    0,0,1,
    0,-1,0;
  Eigen::MatrixXi SFgt(5,3);
  SFgt<<
  0,3,4,
  0,4,5,
  0,5,3,
  1,6,2,
  2,6,7;
  Eigen::VectorXi SVIgt(8);
  SVIgt<<  0,  0,  1,  1,  2,  3,  4,  3;
  test_common::assert_eq(SV,SVgt);
  test_common::assert_eq(SF,SFgt);
  test_common::assert_eq(SVI,SVIgt);
}

TEST_CASE("split_nonmanifold: vertex-kiss", "[igl]")
{
  using namespace igl;
  Eigen::MatrixXd V(7,3);
  V<<
    0,0,0,
    1,0,0,
    0,1,0,
    0,0,1,
    0,0,2,
    1,0,2,
    0,1,2;
  Eigen::MatrixXi F(6,3);
  F<<
    0,1,3,
    1,2,3,
    2,0,3,
    4,5,3,
    5,6,3,
    6,4,3;
  Eigen::MatrixXd SV;
  Eigen::MatrixXi SF;
  Eigen::VectorXi SVI;
  igl::split_nonmanifold(V,F,SV,SF,SVI);
  Eigen::MatrixXd SVgt(8,3);
  SVgt<<
    0,0,0,
    1,0,0,
    0,1,0,
    0,0,2,
    1,0,2,
    0,1,2,
    0,0,1,
    0,0,1;
  Eigen::MatrixXi SFgt(6,3);
  SFgt<<
    0,1,6,
    1,2,6,
    2,0,6,
    3,4,7,
    4,5,7,
    5,3,7;
  Eigen::VectorXi SVIgt(8);
  SVIgt<<  0,  1,  2,  4,  5,  6,  3,  3;
  test_common::assert_eq(SV,SVgt);
  test_common::assert_eq(SF,SFgt);
  test_common::assert_eq(SVI,SVIgt);
}

TEST_CASE("split_nonmanifold: non-orientable", "[igl]")
{
  using namespace igl;
  Eigen::MatrixXd V(6,3);
  V<<
     6, 0, 0,
     4, 0, 0,
    -3, 5, 0,
    -2, 4, 0,
    -2,-4, 1,
    -3,-5,-1;
  Eigen::MatrixXi F(6,3);
  F<<
    0,2,1,
    2,3,1,
    2,4,3,
    4,5,3,
    4,1,5,
    1,0,5;
  Eigen::MatrixXd SV;
  Eigen::MatrixXi SF;
  Eigen::VectorXi SVI;
  igl::split_nonmanifold(V,F,SV,SF,SVI);
  Eigen::MatrixXd SVgt(8,3);
  SVgt<<
    6,0,0,
    -3,5,0,
    -2,-4,1,
    4,0,0,
    -2,4,0,
    -3,-5,-1,
    6,0,0,
    4,0,0;
  Eigen::MatrixXi SFgt(6,3);
  SFgt<<
    0,1,7,
    1,4,7,
    1,2,4,
    2,5,4,
    2,3,5,
    3,6,5;
  Eigen::VectorXi SVIgt(8);
  SVIgt<<  0,  2,  4,  1,  3,  5,  0,  1;
  test_common::assert_eq(SV,SVgt);
  test_common::assert_eq(SF,SFgt);
  test_common::assert_eq(SVI,SVIgt);
}

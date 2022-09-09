#include <test_common.h>
#include <igl/signed_distance.h>

TEST_CASE("signed_distance: single_tet", "[igl]")
{
  Eigen::MatrixXd V(4,3);
  V<<
    0,0,0,
    1,0,0,
    0,1,0,
    0,0,1;
  Eigen::MatrixXi F(4,3);
  F<<
    0,1,3,
    0,2,1,
    0,3,2,
    1,2,3;

  Eigen::MatrixXd P(1,3);
  P<<0.5,0.5,0.5;
  for(const igl::SignedDistanceType type :
      {
      igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL  ,
      igl::SIGNED_DISTANCE_TYPE_WINDING_NUMBER,
      igl::SIGNED_DISTANCE_TYPE_DEFAULT       ,
      igl::SIGNED_DISTANCE_TYPE_UNSIGNED      ,
      igl::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER 
      })
  {
    Eigen::VectorXd S;
    Eigen::VectorXi I;
    Eigen::MatrixXd C,N;
    igl::signed_distance( P,V,F,type,S,I,C,N);
    Eigen::VectorXd Sexact (1,1);Sexact<<sqrt(1./12.);

    if (type == igl::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER) {
      // loosen tolerance on fwn. 
      test_common::assert_near(S,Sexact,1e-7);
    } else {
      test_common::assert_near(S,Sexact,1e-15);
    }
    
  }
}

TEST_CASE("signed_distance: single_triangle", "[igl]")
{
  Eigen::MatrixXd V(3,2);
  V<<
    0,0,
    1,0,
    0,1;
  Eigen::MatrixXi F(3,2);
  F<<
    0,1,
    1,2,
    2,0;
  Eigen::MatrixXd P(1,2);
  P<<1,1;
  for(const igl::SignedDistanceType type :
      {
      igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL  ,
      igl::SIGNED_DISTANCE_TYPE_WINDING_NUMBER,
      igl::SIGNED_DISTANCE_TYPE_DEFAULT       ,
      igl::SIGNED_DISTANCE_TYPE_UNSIGNED
      })
  {
    Eigen::VectorXd S;
    Eigen::VectorXi I;
    Eigen::MatrixXd C,N;
    igl::signed_distance( P,V,F,type,S,I,C,N);
    Eigen::VectorXd Sexact (1,1);Sexact<<sqrt(2.)/2.;
    test_common::assert_near(S,Sexact,1e-15);
  }
}


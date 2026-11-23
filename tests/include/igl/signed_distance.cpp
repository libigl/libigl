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

TEST_CASE("signed_distance: dimension-templates", "[igl]")
{
  Eigen::MatrixXd VX3(3,3);
  VX3<<
    0,0,0,
    1,0,0,
    0,1,0;
  Eigen::MatrixXi F(1,3);
  F<<0,1,2;
  Eigen::MatrixXi E(3,2);
  E<<
    0,1,
    1,2,
    2,0;
  Eigen::MatrixXd VX2 = VX3.leftCols(2);
  Eigen::MatrixX2d V2 = VX3.leftCols(2);
  Eigen::MatrixX3d V3 = VX3.leftCols(3);

  Eigen::MatrixXd PX3 = VX3;
  Eigen::MatrixXd PX2 = VX2;
  Eigen::MatrixX2d P2 = VX3.leftCols(2);
  Eigen::MatrixX3d P3 = VX3.leftCols(3);
  Eigen::VectorXd S;
  Eigen::VectorXi I;
  Eigen::MatrixXd C,N;
  Eigen::MatrixX2d C2,N2;
  Eigen::MatrixX3d C3,N3;
  igl::SignedDistanceType type = igl::SIGNED_DISTANCE_TYPE_DEFAULT;
  double lower_bound = std::numeric_limits<double>::min();
  double upper_bound = std::numeric_limits<double>::max();


  Eigen::MatrixX2i E2 = E;
  Eigen::MatrixX3i F3 = F;
  igl::signed_distance(P2,V2,E2,type,lower_bound,upper_bound,S,I,C2,N2);
  igl::signed_distance(P3,V3,F3,type,lower_bound,upper_bound,S,I,C3,N3);

  igl::signed_distance(PX2,VX2,E,type,lower_bound,upper_bound,S,I,C,N);
  igl::signed_distance(PX3,VX3,F,type,lower_bound,upper_bound,S,I,C,N);
}

TEMPLATE_TEST_CASE("signed_distance: pseudonormal scalar", "[igl]", float, double)
{
  // Exercise closest faces, edges, and vertices in both area branches.
  using Scalar = TestType;
  using Matrix = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>;
  Eigen::MatrixXi F(4,3);
  F << 0,1,3, 0,2,1, 0,3,2, 1,2,3;
  for(const Scalar scale : {Scalar(1),Scalar(1)/Scalar(1024)})
  {
    Matrix V(4,3);
    V << 0,0,0, 1,0,0, 0,1,0, 0,0,1;
    V *= scale;
    Matrix P(4,3);
    P << Scalar(0.125),Scalar(0.125),Scalar(0.125),
         Scalar(0.25),Scalar(0.25),-1,
         Scalar(0.5),-1,-1,
         -1,-1,-1;
    P *= scale;
    Eigen::Matrix<Scalar,Eigen::Dynamic,1> S;
    Eigen::VectorXi I;
    Matrix C,N;
    igl::signed_distance(P,V,F,igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL,S,I,C,N);
    Eigen::Matrix<Scalar,4,1> expected;
    expected << -Scalar(0.125),Scalar(1),std::sqrt(Scalar(2)),std::sqrt(Scalar(3));
    expected *= scale;
    test_common::assert_near(S,expected,Scalar(16)*std::numeric_limits<Scalar>::epsilon()*scale);
    REQUIRE(C.allFinite());
    REQUIRE(N.allFinite());

    // Preserve typed bounds and out-of-range NaNs.
    igl::signed_distance(P,V,F,igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL,
      -scale/Scalar(4),scale/Scalar(4),S,I,C,N);
    REQUIRE(S(0) == expected(0));
    REQUIRE(std::isnan(S(1)));
  }
}

TEMPLATE_TEST_CASE("signed_distance: 2d pseudonormal scalar", "[igl]", float, double)
{
  using Scalar = TestType;
  using Matrix = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>;
  Matrix V(3,2);
  V << 0,0, 1,0, 0,1;
  Eigen::MatrixXi E(3,2);
  E << 0,1, 1,2, 2,0;
  Matrix P(3,2);
  P << Scalar(0.25),Scalar(0.25), Scalar(0.5),-1, -1,-1;
  Eigen::Matrix<Scalar,Eigen::Dynamic,1> S;
  Eigen::VectorXi I;
  Matrix C,N;
  igl::signed_distance(P,V,E,igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL,S,I,C,N);
  Eigen::Matrix<Scalar,3,1> expected;
  expected << -Scalar(0.25),Scalar(1),std::sqrt(Scalar(2));
  test_common::assert_near(S,expected,Scalar(16)*std::numeric_limits<Scalar>::epsilon());
  REQUIRE(C.allFinite());
  REQUIRE(N.allFinite());
}

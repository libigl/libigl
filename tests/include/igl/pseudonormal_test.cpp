// Exercise scalar/shape combinations beyond the static library instantiations.
#ifdef IGL_STATIC_LIBRARY
#  undef IGL_STATIC_LIBRARY
#endif

#include <test_common.h>
#include <igl/pseudonormal_test.h>
#include <igl/signed_distance.h>

TEMPLATE_TEST_CASE("pseudonormal_test: scalar feature selection", "[igl]", float, double)
{
  using Scalar = TestType;
  using Row = Eigen::Matrix<Scalar,1,3>;
  using Matrix = Eigen::Matrix<Scalar,Eigen::Dynamic,3>;
  Eigen::Matrix<int,1,3> F;
  F << 0,1,2;
  Eigen::Vector3i EMAP;
  EMAP << 0,1,2;
  Matrix FN(1,3),VN(3,3),EN(3,3);
  FN << 0,0,1;
  VN << 1,0,0, 1,0,0, 1,0,0;
  EN << 0,1,0, 0,1,0, 0,1,0;
  for(const Scalar scale : {Scalar(1),Scalar(1)/Scalar(1024)})
  {
    Matrix V(3,3);
    V << 0,0,0, 1,0,0, 0,1,0;
    V *= scale;
    Matrix closest(3,3);
    closest << 0,0,0, Scalar(0.5),0,0, Scalar(0.25),Scalar(0.25),0;
    closest *= scale;
    for(int feature = 0; feature < 3; ++feature)
    {
      Row expected = Row::Zero();
      expected(feature) = 1;
      for(const Scalar direction : {Scalar(-1),Scalar(1)})
      {
        Row c = closest.row(feature);
        const Row q = c + direction*expected;
        Row n;
        Scalar s;
        igl::pseudonormal_test(V,F,FN,VN,EN,EMAP,q,0,c,s,n);
        test_common::assert_eq(n,expected);
        REQUIRE(s == direction);
      }
    }
  }
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

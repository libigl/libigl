#include <test_common.h>
#include <igl/pseudonormal_test.h>

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

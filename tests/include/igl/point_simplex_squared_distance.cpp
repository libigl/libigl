// Exercise scalar/shape combinations beyond the static library instantiations.
#ifdef IGL_STATIC_LIBRARY
#  undef IGL_STATIC_LIBRARY
#endif

#include <test_common.h>
#include <igl/point_simplex_squared_distance.h>

TEST_CASE("point_simplex_squared_distance: float barycentric arithmetic", "[igl]")
{
  // A double literal in 1-v-w changes the rounding of the first coordinate.
  Eigen::Matrix<float,3,3> V;
  V << 0,0,0, 1,0,0, 0,1,0;
  Eigen::Matrix<int,1,3> F;
  F << 0,1,2;
  Eigen::RowVector3f p,c,bary;
  p << 1.0f/257.0f,109.0f/257.0f,1.0f;
  float squared_distance;
  igl::point_simplex_squared_distance<3>(p,V,F,0,squared_distance,c,bary);
  const float first_subtraction = 1.0f-bary(1);
  REQUIRE(bary(0) == first_subtraction-bary(2));
  REQUIRE(squared_distance == 1.0f);
}

#include <test_common.h>
#include <igl/nextafter.h>
#include <limits>

TEST_CASE("nextafter: dense", "[igl]" )
{
  Eigen::RowVector4f v, ref;
  v << 0.f, std::numeric_limits<float>::max(), std::numeric_limits<float>::infinity(), std::numeric_limits<float>::quiet_NaN();
  ref << std::numeric_limits<float>::denorm_min(), std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity(), std::numeric_limits<float>::quiet_NaN();

  igl::nextafter(v, 1);

  REQUIRE(v.head(3) == ref.head(3));
  REQUIRE((std::isnan(v(3)) && std::isnan(ref(3))));
}

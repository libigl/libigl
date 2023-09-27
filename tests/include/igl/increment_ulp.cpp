#include <test_common.h>
#include <igl/increment_ulp.h>
#include <limits>

TEST_CASE("increment_ulp: dense", "[igl]" )
{
  Eigen::RowVectorXf v(5), ref(5);
  v << std::numeric_limits<float>::lowest(), 0.f, std::numeric_limits<float>::max(), std::numeric_limits<float>::infinity(), std::numeric_limits<float>::quiet_NaN();
  ref << -std::numeric_limits<float>::infinity(), std::numeric_limits<float>::denorm_min(), std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity(), std::numeric_limits<float>::quiet_NaN();

  igl::increment_ulp(v, 1);

  REQUIRE(v.head(4) == ref.head(4));
  REQUIRE((std::isnan(v(4)) && std::isnan(ref(4))));
}

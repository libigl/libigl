#include <test_common.h>
#include <igl/sign.h>


template <typename T>
void test()
{
  REQUIRE (igl::sign(T(  0)) == 0);
  REQUIRE (igl::sign(T( -0)) == 0);
  REQUIRE (igl::sign(T(  1)) == 1);
  REQUIRE (igl::sign(T( -1)) == -1);
  if constexpr (std::is_floating_point<T>::value)
  {
    REQUIRE (igl::sign(T(  0.0)) == 0);
    REQUIRE (igl::sign(T( -0.0)) == 0);
    REQUIRE (igl::sign(T(0.5)) == 1);
    REQUIRE (igl::sign(T(-0.5)) == -1);
    REQUIRE (igl::sign(T(0.00000000000001)) == 1);
    REQUIRE (igl::sign(T(-0.00000000000001)) == -1);
  }
}

TEST_CASE("sign: cases", "[igl]" )
{
  test<float>();
  test<double>();
  test<int>();
}

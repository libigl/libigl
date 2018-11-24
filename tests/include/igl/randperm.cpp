#include <test_common.h>
#include <igl/randperm.h>
#include <random>

TEST(randperm, default_rng_reproduce_identity)
{
  int n = 100;
  Eigen::VectorXi I1, I2;

  std::srand(6);
  igl::randperm(100, I1);
  std::srand(6);
  igl::randperm(100, I2);

  test_common::assert_eq(I1, I2);
}

TEST(randperm, custom_rng_reproduce_identity)
{
  int n = 100;
  Eigen::VectorXi I1, I2;
  std::minstd_rand rng1(6);
  std::minstd_rand rng2(6);

  igl::randperm(100, I1, rng1.min(), rng1.max(),
                [&rng1]()->int64_t { return rng1(); });
  igl::randperm(100, I2, rng2.min(), rng2.max(),
                [&rng2]()->int64_t { return rng2(); });

  test_common::assert_eq(I1, I2);
}

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

namespace randperm
{
  template<typename URBG>
  void test_reproduce()
  {
    int n = 100;
    Eigen::VectorXi I1, I2;
    URBG rng1(6);
    URBG rng2(6);

    igl::randperm(100, I1, rng1);
    igl::randperm(100, I2, rng2);

    test_common::assert_eq(I1, I2);
  }
}

TEST(randperm, minstd_rand0_reproduce_identity)
{
  randperm::test_reproduce<std::minstd_rand0>();
}
TEST(randperm, minstd_rand_reproduce_identity)
{
  randperm::test_reproduce<std::minstd_rand>();
}
TEST(randperm, mt19937_reproduce_identity)
{
  randperm::test_reproduce<std::mt19937>();
}
TEST(randperm, mt19937_64_reproduce_identity)
{
  randperm::test_reproduce<std::mt19937_64>();
}
TEST(randperm, ranlux24_base_reproduce_identity)
{
  randperm::test_reproduce<std::ranlux24_base>();
}
TEST(randperm, ranlux48_base_reproduce_identity)
{
  randperm::test_reproduce<std::ranlux48_base>();
}
TEST(randperm, ranlux24_reproduce_identity)
{
  randperm::test_reproduce<std::ranlux24>();
}
TEST(randperm, ranlux48_reproduce_identity)
{
  randperm::test_reproduce<std::ranlux48>();
}
TEST(randperm, knuth_b_reproduce_identity)
{
  randperm::test_reproduce<std::knuth_b>();
}

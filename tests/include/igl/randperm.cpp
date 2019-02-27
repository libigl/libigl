#include <test_common.h>
#include <igl/randperm.h>
#include <random>

TEST_CASE("randperm: default_rng_reproduce_identity", "[igl]")
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
    Eigen::MatrixXi Ix1, Ix2;
    URBG rng1(6);
    URBG rng2(6);

    igl::randperm(100, I1, rng1);
    igl::randperm(100, I2, rng2);

    igl::randperm(100, Ix1, rng1);
    igl::randperm(100, Ix2, rng2);

    test_common::assert_eq(I1, I2);
    test_common::assert_eq(Ix1, Ix2);

    test_common::assert_neq(I1, Ix1);
    test_common::assert_neq(I2, Ix2);
  }
}

TEST_CASE("randperm: minstd_rand0_reproduce_identity", "[igl]")
{
  randperm::test_reproduce<std::minstd_rand0>();
}
TEST_CASE("randperm: minstd_rand_reproduce_identity", "[igl]")
{
  randperm::test_reproduce<std::minstd_rand>();
}
TEST_CASE("randperm: mt19937_reproduce_identity", "[igl]")
{
  randperm::test_reproduce<std::mt19937>();
}
TEST_CASE("randperm: mt19937_64_reproduce_identity", "[igl]")
{
  randperm::test_reproduce<std::mt19937_64>();
}
TEST_CASE("randperm: ranlux24_base_reproduce_identity", "[igl]")
{
  randperm::test_reproduce<std::ranlux24_base>();
}
TEST_CASE("randperm: ranlux48_base_reproduce_identity", "[igl]")
{
  randperm::test_reproduce<std::ranlux48_base>();
}
TEST_CASE("randperm: ranlux24_reproduce_identity", "[igl]")
{
  randperm::test_reproduce<std::ranlux24>();
}
TEST_CASE("randperm: ranlux48_reproduce_identity", "[igl]")
{
  randperm::test_reproduce<std::ranlux48>();
}
TEST_CASE("randperm: knuth_b_reproduce_identity", "[igl]")
{
  randperm::test_reproduce<std::knuth_b>();
}
TEST_CASE("randperm: default_identity", "[igl]")
{
  int n = 100;
  Eigen::VectorXi I1, I2;
  Eigen::MatrixXi Ix1, Ix2;

  std::srand(0);
  igl::randperm(100, I1);
  igl::randperm(100, Ix1);
  std::srand(0);
  igl::randperm(100, I2);
  igl::randperm(100, Ix2);

  test_common::assert_eq(I1, I2);
  test_common::assert_eq(Ix1, Ix2);
}

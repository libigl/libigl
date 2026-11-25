#include <test_common.h>
#include <igl/median.h>
#include <algorithm>
#include <random>
#include <vector>

namespace
{
  // Median computed by fully sorting a copy of the input.
  double sorted_median(const Eigen::VectorXd & Z)
  {
    std::vector<double> sorted(Z.data(),Z.data()+Z.size());
    std::sort(sorted.begin(),sorted.end());
    const size_t n = sorted.size()/2;
    return sorted.size()%2 == 0 ? 0.5*(sorted[n]+sorted[n-1]) : sorted[n];
  }
}

TEST_CASE("median: odd_and_even", "[igl]")
{
  double m;
  Eigen::VectorXd odd(3);
  odd<<3,1,2;
  REQUIRE(igl::median(odd,m));
  REQUIRE(m == 2);
  Eigen::VectorXd even(4);
  even<<4,1,3,2;
  REQUIRE(igl::median(even,m));
  REQUIRE(m == 2.5);
  REQUIRE(!igl::median(Eigen::VectorXd(),m));
}

TEST_CASE("median: tied_keys", "[igl]")
{
  // Sorted: 0 0 0 0 0 0 1 2 2 3. Both middle values are 0.
  Eigen::VectorXd Z(10);
  Z<<2,3,1,0,0,0,2,0,0,0;
  double m;
  REQUIRE(igl::median(Z,m));
  REQUIRE(m == 0);
}

TEST_CASE("median: matches_sorted_reference", "[igl]")
{
  std::mt19937 rng(0);
  for(const int n : {1,2,3,10,11,128})
  {
    // Few distinct values produce many ties; many produce almost none.
    for(const unsigned int num_values : {4u,1000000u})
    {
      for(int trial = 0; trial < 2000; trial++)
      {
        Eigen::VectorXd Z(n);
        for(int i = 0; i < n; i++) { Z(i) = double(rng()%num_values); }
        double m;
        REQUIRE(igl::median(Z,m));
        REQUIRE(m == sorted_median(Z));
      }
    }
  }
}

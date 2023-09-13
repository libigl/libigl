#include <test_common.h>
#include <igl/adjacency_list.h>

TEST_CASE("adjacency_list: simple", "[igl]")
{
  for(int off = 0;off<2;off++)
  {
    Eigen::MatrixXi F(2,3);
    F << 0,1,2,
         0,2,3;
    F.array() += off;
    std::vector<std::vector<int> > A;
    igl::adjacency_list(F,A,true);
    REQUIRE(A.size() == 4+off);
    REQUIRE(A[0+off].size() == 3);
    REQUIRE(A[1+off].size() == 2);
    REQUIRE(A[2+off].size() == 3);
    REQUIRE(A[3+off].size() == 2);
    REQUIRE(A[0+off][0]     == 1+off);
    REQUIRE(A[0+off][1]     == 2+off);
    REQUIRE(A[0+off][2]     == 3+off);
    REQUIRE(A[1+off][0]     == 2+off);
    REQUIRE(A[1+off][1]     == 0+off);
    REQUIRE(A[2+off][0]     == 3+off);
    REQUIRE(A[2+off][1]     == 0+off);
    REQUIRE(A[2+off][2]     == 1+off);
    REQUIRE(A[3+off][0]     == 0+off);
    REQUIRE(A[3+off][1]     == 2+off);
  }
}

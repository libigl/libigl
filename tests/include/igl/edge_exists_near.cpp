#include <test_common.h>
#include <igl/edge_exists_near.h>
#include <igl/unique_edge_map.h>

TEST_CASE("edge_exists_near: tet", "[igl]")
{
  const Eigen::MatrixXi F = (Eigen::MatrixXi(4,3)<<
     0,1,2,
     0,2,3,
     0,3,1,
     1,3,2).finished();
  Eigen::MatrixXi E,uE;
  Eigen::VectorXi EMAP;
  std::vector<std::vector<int> > uE2E;
  igl::unique_edge_map(F,E,uE,EMAP,uE2E);
  for(int uei = 0;uei<uE2E.size();uei++)
  {
    for(int i = 0;i<4;i++)
    {
      for(int j = 0;j<4;j++)
      {
        if(i != j)
        {
          REQUIRE (igl::edge_exists_near(uE,EMAP,uE2E,i,j,uei));
        }
      }
    }
  }
}

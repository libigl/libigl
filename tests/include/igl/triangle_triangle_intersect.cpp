#include "test_common.h"
#include <igl/triangle_triangle_intersect.h>
#include <igl/edge_flaps.h>

TEST_CASE("triangle_triangle_intersect: shared-edge", "[igl]" )
{
  Eigen::MatrixXd V(4,3);
  V<<
    0,0,0,
    1,0,0,
    0,1,0,
    -1,0,0;
  Eigen::MatrixXi F(2,3);
  F <<
    0,1,2,
    0,2,3;
  //      2
  //     /|\
  //    / | \
  //   /  |  \
  //  3---0---1
  const int f = 0;
  const int c = 1;
  const int g = 1;

  Eigen::MatrixXi E,EF,EI;
  Eigen::VectorXi EMAP;
  igl::edge_flaps(F,E,EMAP,EF,EI);


  bool ret;
  ret = igl::triangle_triangle_intersect(V,F,E,EMAP,EF,f,c,V.row(F(f,c)),g);
  REQUIRE(ret == false);

  for(const double epsilon : {0.,1e-15,-1e-15})
  {
    //  2
    //  |\⟍
    //  | \ ⟍
    //  |  \  ⟍
    //  0---1---3
    V.row(3) << 2.0+epsilon,0,0;
    ret = igl::triangle_triangle_intersect(V,F,E,EMAP,EF,f,c,V.row(F(f,c)),g);
    REQUIRE(ret == true);
  }

}

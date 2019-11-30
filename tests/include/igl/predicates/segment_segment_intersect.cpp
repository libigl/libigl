#include <test_common.h>
#include <igl/predicates/segment_segment_intersect.h>
#include <iomanip>

TEST_CASE("segment_segment_intersect: robust", "[igl/predicates]")
{
  // example 1: vanila intersecting case
  auto A1 = Eigen::RowVector2d(0, 128.5);
  auto B1 = Eigen::RowVector2d(-77.44,1.2);
  auto C1 = Eigen::RowVector2d(-83.2,2.8);
  auto D1 = Eigen::RowVector2d(-1.0,-1.0);

  bool check1 = igl::predicates::segment_segment_intersect(A1,B1,C1,D1);
  REQUIRE(check1 == true);
  
  // example 2: colinear overlapping
  auto A2 = Eigen::RowVector2d(1.0,5.0);
  auto B2 = Eigen::RowVector2d(1.0,9.0);
  auto C2 = Eigen::RowVector2d(1.0,8.0);
  auto D2 = Eigen::RowVector2d(1.0,12.0);
  
  bool check2 = igl::predicates::segment_segment_intersect(A2,B2,C2,D2);
  REQUIRE(check2 == true);

  // example 3: colinear touching endpoint
  auto A3 = Eigen::RowVector2d(0.0,0.0);
  auto B3 = Eigen::RowVector2d(1.5,1.5);
  auto C3 = Eigen::RowVector2d(1.5,1.5);
  auto D3 = Eigen::RowVector2d(2.0,2.0);
  
  bool check3 = igl::predicates::segment_segment_intersect(A3,B3,C3,D3);
  REQUIRE(check3 == true);

  // example 6: colinear not touching endpoint
  double eps = 1e-14;
  auto A4 = Eigen::RowVector2d(0.0,0.0);
  auto B4 = Eigen::RowVector2d(1.5,1.5);
  auto C4 = Eigen::RowVector2d(1.5+eps,1.5+eps);
  auto D4 = Eigen::RowVector2d(2.0,2.0);
  bool check4 = igl::predicates::segment_segment_intersect(A4,B4,C4,D4);
  REQUIRE(check4 == false);
  
}
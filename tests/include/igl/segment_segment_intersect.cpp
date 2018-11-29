#include <test_common.h>
#include <igl/segment_segment_intersect.h>

TEST_CASE("segment_segment_intersect: examples", "[igl]")
{
  // example 1
  // https://github.com/libigl/libigl/issues/965

  auto s1 = Eigen::RowVector3d(581, 388, -54);
  auto dir1 = Eigen::RowVector3d(0.75, -4.36, 0.64);

  auto s2 = Eigen::RowVector3d(636, 77, -10);
  auto dir2 = Eigen::RowVector3d(-2.79, 0.96, 0.39);

  double a1, a2;
  bool sect1 = igl::segments_intersect(s1, dir1, s2, dir2, a1, a2);

  bool intersectCondition1 = (a1 >= 0 && a1 <= 1 && a2 >= 0 && a2 <= 1);

  REQUIRE(sect1 == false);
  REQUIRE(intersectCondition1 == false);

  // example 2
  // https://github.com/libigl/libigl/issues/957

  auto s3 = Eigen::RowVector3d(-56.6, 0, -201.7);
  auto dir3 = Eigen::RowVector3d(0, 0, 1);

  auto s4 = Eigen::RowVector3d(-65.6, 0, 258.8);
  auto dir4 = Eigen::RowVector3d(15.4, 0, 1.9);

  double a3, a4;
  bool sect2 = igl::segments_intersect(s3, dir3, s4, dir4, a3, a4);

  bool intersectCondition2 = (a3 >= 0 && a3 <= 1 && a4 >= 0 && a4 <= 1);

  REQUIRE(sect2 == false);
  REQUIRE(intersectCondition2 == false);
}

#include <test_common.h>
#include <igl/predicates/point_in_convex_hull.h>

TEST_CASE("point_in_convex_hull: simple", "[igl/predicates]")
{
  Eigen::RowVector2d a(0,0);
  Eigen::RowVector2d b(1,1);
  Eigen::RowVector2d c(2,-1);
  Eigen::RowVector2d d(3,0);

  {
    Eigen::RowVector2d q(1.5,0);
    auto r = igl::predicates::point_in_convex_hull(q,a,b,c,d);
    REQUIRE(r == igl::Orientation::POSITIVE);
  }
  {
    Eigen::RowVector2d q(-1,0);
    auto r = igl::predicates::point_in_convex_hull(q,a,b,c,d);
    REQUIRE(r == igl::Orientation::NEGATIVE);
  }
  {
    Eigen::RowVector2d q(0,0);
    auto r = igl::predicates::point_in_convex_hull(q,a,b,c,d);
    REQUIRE(r == igl::Orientation::COLLINEAR);
  }
  {
    Eigen::RowVector2d q(0.5,0.5);
    auto r = igl::predicates::point_in_convex_hull(q,a,b,c,d);
    REQUIRE(r == igl::Orientation::COLLINEAR);
  }
  {
    Eigen::RowVector2d q(1.0,-0.4);
    auto r = igl::predicates::point_in_convex_hull(q,a,b,c,d);
    REQUIRE(r == igl::Orientation::POSITIVE);
  }
  {
    Eigen::RowVector2d q(2.0,0.6);
    auto r = igl::predicates::point_in_convex_hull(q,a,b,c,d);
    REQUIRE(r == igl::Orientation::NEGATIVE);
  }
}

TEST_CASE("point_in_convex_hull: degenerate", "[igl/predicates]")
{
  // Triangle
  {
    Eigen::RowVector2d a(0,0);
    Eigen::RowVector2d b(1,0);
    Eigen::RowVector2d c(0,1);
    Eigen::RowVector2d d(0.1,0.1);
    REQUIRE(igl::predicates::point_in_convex_hull(Eigen::RowVector2d(0.1,0.1),a,b,c,d) == igl::Orientation::POSITIVE);
    REQUIRE(igl::predicates::point_in_convex_hull(Eigen::RowVector2d(0.2,0.2),a,b,c,d) == igl::Orientation::POSITIVE);
    REQUIRE(igl::predicates::point_in_convex_hull(Eigen::RowVector2d(-0.2,-0.2),a,b,c,d) == igl::Orientation::NEGATIVE);
    REQUIRE(igl::predicates::point_in_convex_hull(Eigen::RowVector2d(0.5,0.5),a,b,c,d) == igl::Orientation::COLLINEAR);
  }

  // Segment
  {
    Eigen::RowVector2d a(0,0);
    Eigen::RowVector2d b(1,0);
    Eigen::RowVector2d c(0.5,0);
    Eigen::RowVector2d d(0.25,0);
    REQUIRE(igl::predicates::point_in_convex_hull(Eigen::RowVector2d(0.25,0),a,b,c,d) == igl::Orientation::COLLINEAR);
    REQUIRE(igl::predicates::point_in_convex_hull(Eigen::RowVector2d(0.5,0),a,b,c,d) == igl::Orientation::COLLINEAR);
    REQUIRE(igl::predicates::point_in_convex_hull(Eigen::RowVector2d(1.5,0),a,b,c,d) == igl::Orientation::NEGATIVE);
    REQUIRE(igl::predicates::point_in_convex_hull(Eigen::RowVector2d(-0.5,0),a,b,c,d) == igl::Orientation::NEGATIVE);
    REQUIRE(igl::predicates::point_in_convex_hull(Eigen::RowVector2d(0.5,1),a,b,c,d) == igl::Orientation::NEGATIVE);
  }

  // Point
  {
    Eigen::RowVector2d a(0,0);
    Eigen::RowVector2d b(0,0);
    Eigen::RowVector2d c(0,0);
    Eigen::RowVector2d d(0,0);
    REQUIRE(igl::predicates::point_in_convex_hull(Eigen::RowVector2d(0,0),a,b,c,d) == igl::Orientation::COLLINEAR);
    REQUIRE(igl::predicates::point_in_convex_hull(Eigen::RowVector2d(1,0),a,b,c,d) == igl::Orientation::NEGATIVE);
    REQUIRE(igl::predicates::point_in_convex_hull(Eigen::RowVector2d(0,1),a,b,c,d) == igl::Orientation::NEGATIVE);
    REQUIRE(igl::predicates::point_in_convex_hull(Eigen::RowVector2d(-1,0),a,b,c,d) == igl::Orientation::NEGATIVE);
    REQUIRE(igl::predicates::point_in_convex_hull(Eigen::RowVector2d(0,-1),a,b,c,d) == igl::Orientation::NEGATIVE);
  }

}

TEST_CASE("point_in_convex_hull: double-degenerate", "[igl/predicates]")
{
  Eigen::RowVector2d q(-0.86700610997963345,-0.50040733197556009);
  Eigen::Matrix<double,4,2> C(4,2);
  C<< 
    -0.44999999999999996,0,
    -0.29999999999999999,0.17999999999999999,
    -0.14999999999999999,0.35999999999999999,
    -0,0.54000000000000004;
  const auto res = igl::predicates::point_in_convex_hull(
      q,
      C.row(0).eval(),
      C.row(1).eval(),
      C.row(2).eval(),
      C.row(3).eval());
  REQUIRE( res == igl::Orientation::NEGATIVE );
}


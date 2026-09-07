#include <test_common.h>
#include <igl/copyleft/cgal/snap_rounding.h>

TEST_CASE("igl_copyleft_cgal_snap_rounding: crossing_segments", "[igl/copyleft/cgal]")
{
  // Two segments crossing exactly at the integer point (5,5). Snap rounding
  // should insert (5,5) as a Steiner vertex and split both segments through it.
  // This exercises the "centroid of hits" path where a segment passes through a
  // hot pixel.
  const Eigen::MatrixXd V = (Eigen::MatrixXd(4,2)<<
     0, 0,
    10,10,
     0,10,
    10, 0).finished();
  const Eigen::MatrixXi E = (Eigen::MatrixXi(2,2)<<
    0,1,
    2,3).finished();
  Eigen::MatrixXd VI;
  Eigen::MatrixXi EI,J;
  igl::copyleft::cgal::snap_rounding(V,E,VI,EI,J);

  // Output vertices are integer-valued, and the rounded copies of V come first.
  REQUIRE(VI.cols() == 2);
  REQUIRE(VI.rows() == 5);
  REQUIRE((VI.array() == VI.array().round()).all());
  REQUIRE(VI.topRows(V.rows()) == V);
  // The lone new vertex is the intersection (5,5).
  REQUIRE(VI(4,0) == 5);
  REQUIRE(VI(4,1) == 5);

  // Each of the two input segments is split at the intersection => 4 edges, and
  // every output edge is incident on the inserted vertex 4.
  REQUIRE(EI.cols() == 2);
  REQUIRE(EI.rows() == 4);
  for(int i = 0;i<EI.rows();i++)
  {
    REQUIRE((EI(i,0) == 4 || EI(i,1) == 4));
    REQUIRE(EI(i,0) != EI(i,1));
  }

  // J maps each output edge back to a valid parent segment, two children each.
  REQUIRE(J.rows() == EI.rows());
  REQUIRE((J.array() >= 0).all());
  REQUIRE((J.array() < E.rows()).all());
  REQUIRE((J.array() == 0).count() == 2);
  REQUIRE((J.array() == 1).count() == 2);
}

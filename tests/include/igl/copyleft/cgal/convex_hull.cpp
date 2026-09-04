#include <test_common.h>
#include <igl/copyleft/cgal/convex_hull.h>

TEST_CASE(
  "igl_copyleft_cgal_convex_hull: cube",
  "[igl/copyleft/cgal/]")
{
  Eigen::MatrixXd V(8,3);
  V << 0,0,0,
       1,0,0,
       1,1,0,
       0,1,0,
       0,0,1,
       1,0,1,
       1,1,1,
       0,1,1;
  Eigen::MatrixXi F;
  igl::copyleft::cgal::convex_hull(V,F);
  REQUIRE(F.rows() == 12); // cube has 12 triangles
}

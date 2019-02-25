#include <test_common.h>
#include <igl/copyleft/cgal/delaunay_triangulation.h>
#include <igl/unique_simplices.h>
#include <igl/matlab_format.h>

TEST_CASE("igl_copyleft_cgal_delaunay_triangulation: two_triangles", "[igl/copyleft/cgal]")
{
  const Eigen::MatrixXd V =
    (Eigen::MatrixXd(4,2)<<
     0,10,
     1,0,
     1,20,
     2,10).finished();
  Eigen::MatrixXi F;
  igl::copyleft::cgal::delaunay_triangulation(V,F);
  // Ground truth
  Eigen::MatrixXi Fgt = (Eigen::MatrixXi(2,3)<<0,1,3,0,3,2).finished();
  REQUIRE (2 == F.rows());
  Eigen::MatrixXi Fu;
  Eigen::VectorXi IA,IC;
  igl::unique_simplices(
    (Eigen::MatrixXi(4,3)<<F,Fgt).finished(),
    Fu,IA,IC);
  // Now new faces w.r.t. ground truth
  REQUIRE (2 == Fu.rows());
}

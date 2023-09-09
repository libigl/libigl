#include <igl/copyleft/cgal/half_space_box.h>
#include <test_common.h>

TEST_CASE("half_space_box: simple", "[igl]")
{
  Eigen::MatrixXd V(2,3);
  V<<0,0,0,
    1,1,1;
  Eigen::RowVector4d equ;
  Eigen::Matrix<CGAL::Epeck::FT,8,3> BV;
  Eigen::Matrix<int,12,3> BF;

  equ << 1,-1,0,0;
  igl::copyleft::cgal::half_space_box(equ,V,BV,BF);
  REQUIRE((BV.colwise().maxCoeff() - BV.colwise().minCoeff()).squaredNorm() > 0);

  equ << 1,1,1,0;
  igl::copyleft::cgal::half_space_box(equ,V,BV,BF);
  REQUIRE((BV.colwise().maxCoeff() - BV.colwise().minCoeff()).squaredNorm() > 0);
}

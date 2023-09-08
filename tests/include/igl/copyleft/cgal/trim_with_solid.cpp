#include <test_common.h>
#include <igl/copyleft/cgal/trim_with_solid.h>
#include <igl/matlab_format.h>
#include <iostream>

TEST_CASE("trim_with_solid: triangle-vs-tet", "[igl/copyleft/cgal]")
{
  // Single Tet
  Eigen::MatrixXd VB(4,3);
  VB<<
    0,0,0,
    1,0,0,
    0,1,0,
    0,0,1;
  Eigen::MatrixXi FB(4,3);
  FB<<
    0,1,3,
    0,2,1,
    0,3,2,
    1,2,3;
  Eigen::MatrixXd VA(3,3);
  VA<< 
    0.2,0.2,0.2,
    0.8,0.2,0.2,
    0.2,0.8,0.2;
  Eigen::MatrixXi FA(1,3);
  FA<<0,1,2;
  using namespace igl::copyleft::cgal;
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::Array<bool,Eigen::Dynamic,1> D;
  Eigen::VectorXi J;
  const auto test = [&]()
  {
    REQUIRE(V.rows() == 5);
    REQUIRE(F.rows() == 3);
    REQUIRE(D.size() == F.rows());
    REQUIRE(J.size() == F.rows());
    REQUIRE(D.count() == 2);
    REQUIRE((J.array() == 0).count() == 3);
  };
  igl::copyleft::cgal::trim_with_solid(VA,FA,VB,FB,CHECK_EACH_FACE,V,F,D,J);
  test();
  igl::copyleft::cgal::trim_with_solid(VA,FA,VB,FB,CHECK_EACH_PATCH,V,F,D,J);
  test();
  igl::copyleft::cgal::trim_with_solid(VA,FA,VB,FB,RESOLVE_BOTH_AND_RESTORE_THEN_CHECK_EACH_PATCH,V,F,D,J);
  test();
}

TEST_CASE("trim_with_solid: two-vs-cube", "[igl/copyleft/cgal]")
{
  Eigen::MatrixXd VA(6,3);
  VA<< 
    3.5301513671875,77.03765869140625,1.3964500427246094,
    3.5709075927734375,77.055023193359375,1.396728515625,
    4.4851531982421875,78.03460693359375,1.6064491271972656,
    3.6009368896484375,77.064865112304688,1.50128173828125,
    3.613006591796875,77.075912475585938,1.5090217590332031,
    4.2685699462890625,77.5335693359375,1.6356697082519531;
  Eigen::MatrixXi FA(2,3);
  FA<< 
    1,2,3,
    0,4,5;
  Eigen::MatrixXd VB(8,3);
  VB<< 
    -1.8674041748046903,-5.6647888183593835,-0.071252098083496196,
    -1.8674041748046903,118.96056518554688,-0.071252098083496196,
    39.21548767089844,-5.6647888183593835,-0.071252098083496196,
    39.21548767089844,118.96056518554688,-0.071252098083496196,
    39.21548767089844,-5.6647888183593835,1.496294059753418,
    39.21548767089844,118.96056518554688,1.496294059753418,
    -1.8674041748046903,-5.6647888183593835,1.496294059753418,
    -1.8674041748046903,118.96056518554688,1.496294059753418;
  Eigen::MatrixXi FB(12,3);
  FB<< 
    1,2,0,
    1,3,2,
    3,4,2,
    3,5,4,
    0,4,6,
    0,2,4,
    7,3,1,
    7,5,3,
    7,0,6,
    7,1,0,
    5,6,4,
    5,7,6;
  using namespace igl::copyleft::cgal;
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::Array<bool,Eigen::Dynamic,1> D;
  Eigen::VectorXi J;
  const auto test = [&]()
  {
    REQUIRE(V.rows() >= VA.rows());
    REQUIRE(FA.rows() == FA.rows());
    REQUIRE(D.size() == F.rows());
    REQUIRE(J.size() == F.rows());
    // some in some out
    REQUIRE( (D.array() == false).count() != 0);
    REQUIRE( (D.array() == true).count() != 0);
  };
  igl::copyleft::cgal::trim_with_solid(VA,FA,VB,FB,CHECK_EACH_FACE,V,F,D,J);
  test();
  igl::copyleft::cgal::trim_with_solid(VA,FA,VB,FB,CHECK_EACH_PATCH,V,F,D,J);
  test();
  igl::copyleft::cgal::trim_with_solid(VA,FA,VB,FB,RESOLVE_BOTH_AND_RESTORE_THEN_CHECK_EACH_PATCH,V,F,D,J);
  test();

}

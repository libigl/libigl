#include <test_common.h>
#include <igl/facet_components.h>

TEST_CASE("facet_components: two_triangles", "[igl]")
{
  const Eigen::MatrixXi F1 = 
    (Eigen::MatrixXi(2,3)<<
     0,1,3,
     0,3,2).finished();
  Eigen::VectorXi C1;
  igl::facet_components(F1,C1);
  REQUIRE(C1(0) == 0);
  REQUIRE(C1(0) == C1(1));

  const Eigen::MatrixXi F2 = 
    (Eigen::MatrixXi(2,3)<<
     0,1,2,
     3,4,5).finished();
  Eigen::VectorXi C2;
  igl::facet_components(F2,C2);
  REQUIRE(C2(0) != C2(1));
  REQUIRE(C2.minCoeff() == 0);
  REQUIRE(C2.maxCoeff() == 1);
}

TEST_CASE("facet_components: truck", "[igl]")
{
  Eigen::MatrixXi F;
  {
    Eigen::MatrixXd V;
    igl::read_triangle_mesh(test_common::data_path("truck.obj"), V, F);
  }
  Eigen::VectorXi C;
  igl::facet_components(F,C);
  REQUIRE(C.maxCoeff()+1 == 59);
}

#include <test_common.h>
#include <igl/principal_curvature.h>
#include <igl/cylinder.h>

TEST_CASE("principal_curvature: cylinder", "[igl]")
{
  using namespace igl;
  const int axis_devisions = 20;
  const int height_devisions = 20;
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  cylinder(axis_devisions,height_devisions,V,F);
  Eigen::MatrixXd PD1,PD2;
  Eigen::VectorXd PV1,PV2;
  //PV1: maximal curvature value for each vertex.
  //PV2: minimal curvature value for each vertex.
  igl::principal_curvature(V,F,PD1,PD2,PV1,PV2);
  REQUIRE (PD1.rows() == V.rows());
  REQUIRE (PD2.rows() == V.rows()); 
  REQUIRE (PD1.cols() == 3);
  REQUIRE (PD2.cols() == 3);
  REQUIRE (PV1.size() == V.rows());
  REQUIRE (PV2.size() == V.rows());
  for(int i = 0; i<PV1.size(); ++i)
  {
    //max curvature is greater than or equal to min curvature
    REQUIRE (PV1[i]>=PV2[i]);
  }
}
#include <test_common.h>
#include <igl/qslim.h>
#include <igl/cylinder.h>
#include <igl/upsample.h>
#include <igl/point_mesh_squared_distance.h>
//#include <igl/hausdorff.h>
#include <igl/writePLY.h>

TEST(qslim,cylinder)
{
  using namespace igl;
  const int axis_devisions = 5;
  const int height_devisions = 2+10;
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  cylinder(axis_devisions,height_devisions,V,F);
  Eigen::MatrixXd U;
  Eigen::MatrixXi G;
  Eigen::VectorXi I,J;
  qslim(V,F,2*axis_devisions,U,G,I,J);
  ASSERT_EQ(axis_devisions*2,U.rows());
  double l,u;
  igl::writePLY("qslim-cylinder-vf.ply",V,F);
  igl::writePLY("qslim-cylinder-ug.ply",U,G);
  const auto & hausdorff_lower_bound = [](
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    Eigen::MatrixXd & U,
    Eigen::MatrixXi & G)->double
  {
    Eigen::MatrixXd Vk;
    Eigen::MatrixXi Fk;
    igl::upsample(V,F,Vk,Fk,5);
    Eigen::MatrixXd C;
    Eigen::VectorXi I;
    Eigen::VectorXd D;
    igl::point_mesh_squared_distance(Vk,U,G,D,I,C);
    return D.array().sqrt().maxCoeff();
  };
  //igl::hausdorff(V,F,U,G,1e-14,l,u);
  ASSERT_NEAR(hausdorff_lower_bound(V,F,U,G),0,2e-10);
  //igl::hausdorff(U,G,V,F,1e-14,l,u);
  ASSERT_NEAR(hausdorff_lower_bound(U,G,V,F),0,2e-10);
}

#include <test_common.h>
#include <igl/internal_angles.h>
#include <igl/internal_angles_intrinsic.h>
#include <igl/matlab_format.h>
#include <igl/squared_edge_lengths.h>

TEST_CASE("internal_angles: 1e-7", "[igl]")
{
  Eigen::MatrixXd V = 
    (Eigen::MatrixXd(3,3)<< 0,0,0,1,1,0,1+1e-7,1,0).finished();

  Eigen::MatrixXi F = (Eigen::MatrixXi(1,3)<<0,1,2).finished();
  Eigen::MatrixXd A;
  igl::internal_angles(V,F,A);

  Eigen::MatrixXd L_sq;
  igl::squared_edge_lengths(V,F,L_sq);
  Eigen::MatrixXd iA;
  igl::internal_angles_intrinsic(L_sq,iA);
  test_common::assert_near(A,iA,1e-6);

  Eigen::MatrixXf fV = V.cast<float>();
  Eigen::MatrixXf fA;
  igl::internal_angles(fV,F,fA);
  test_common::assert_near(A,fA.cast<double>(),1e-6);
  // https://github.com/libigl/libigl/issues/1463

  Eigen::MatrixXf fL_sq;
  igl::squared_edge_lengths(fV,F,fL_sq);
  Eigen::MatrixXf fiA;
  igl::internal_angles_intrinsic(fL_sq,fiA);
  test_common::assert_near(A,fiA.cast<double>(),1e-6);
}

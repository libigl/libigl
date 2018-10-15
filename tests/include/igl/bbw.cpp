#include <gtest/gtest.h>
#include <test_common.h>
#include <igl/boundary_conditions.h>
#include <igl/readMESH.h>
#include <igl/writeDMAT.h>
#include <igl/readTGF.h>
#include <igl/bbw.h>

TEST(bbw, decimated_knight)
{
  Eigen::MatrixXd V,C;
  Eigen::MatrixXi T,F,E;
  igl::readMESH(test_common::data_path("decimated-knight.mesh"),V,T,F);
  igl::readTGF(test_common::data_path("decimated-knight.tgf"),C,E);
  Eigen::MatrixXd W_groundtruth, Was, Wmo;
  igl::readDMAT(
    test_common::data_path("decimated-knight-matlab-active-set.dmat"),W_groundtruth);
  Eigen::VectorXi b;
  Eigen::MatrixXd bc;
  igl::boundary_conditions(V,T,C,Eigen::VectorXi(),E,Eigen::MatrixXi(),b,bc);
  igl::BBWData params;
  params.active_set_params.max_iter = 100;
  igl::bbw(V,T,b,bc,params,Was);
  // igl::writeDMAT("decimated-knight-as.dmat",Was);
  ASSERT_LT( (Was-W_groundtruth).array().abs().maxCoeff() ,1e-4);
}


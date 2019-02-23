#include <test_common.h>
#include <igl/boundary_conditions.h>
#include <igl/readMESH.h>
#include <igl/writeDMAT.h>
#include <igl/readTGF.h>
#include <igl/mosek/bbw.h>

TEST_CASE("mosek_bbw: decimated_knight", "[igl/copyleft/mosek]")
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
  igl::mosek::MosekData mosek_params;
  igl::mosek::bbw(V,T,b,bc,params,mosek_params,Wmo);
  igl::writeDMAT("decimated-knight-mo.dmat",Wmo);
  // Mosek is less accurate
  REQUIRE (1e-3 > (Wmo-W_groundtruth).array().abs().maxCoeff());
}

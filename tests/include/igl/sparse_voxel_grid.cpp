#include <test_common.h>
#include <igl/sparse_voxel_grid.h>
#include <igl/unique_rows.h>

TEST_CASE("sparse_voxel_grid: unique", "[igl]" )
{
  const std::function<double(const Eigen::RowVector3d & x)> f = 
    [&](const Eigen::RowVector3d & x)->double
  {
    return x.norm() - 1.0;
  };
  Eigen::RowVector3d p0(0,1.0,0);
  Eigen::MatrixXd GV;
  Eigen::VectorXd Gf;
  Eigen::Matrix<int,Eigen::Dynamic,8> GI;
  igl::sparse_voxel_grid(p0,f,1,1024,Gf,GV,GI);
  Eigen::MatrixXd uGV;
  Eigen::VectorXi _1,_2;
  igl::unique_rows(GV,uGV,_1,_2);
  REQUIRE(GV.rows() == uGV.rows());
}


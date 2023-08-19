
#include <Eigen/Core>

namespace igl
{
  class ARAPData;
  template <typename LbsMatrixType, typename SSCALAR>
  // Seriously, windows? Why are struct and class different?
  struct ArapDOFData;
}

void precomputation(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & W,
  Eigen::MatrixXd & M,
  Eigen::VectorXi & b,
  Eigen::MatrixXd & L,
  igl::ARAPData & arap_data,
  igl::ARAPData & arap_grouped_data,
  igl::ArapDOFData<Eigen::MatrixXd,double> & arap_dof_data);

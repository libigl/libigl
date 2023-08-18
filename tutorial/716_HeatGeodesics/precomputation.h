#include <Eigen/Core>

namespace igl
{
  template <typename Scalar>
  class HeatGeodesicsData;
}

bool precomputation(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const double t,
  igl::HeatGeodesicsData<double>& data,
    );

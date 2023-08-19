#include <Eigen/Core>
namespace igl
{
  template <typename Scalar>
  struct HeatGeodesicsData;
}

bool update(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const double t,
  const double x,
  const double y,
  const Eigen::Matrix4f& model,
  const Eigen::Matrix4f& proj,
  const Eigen::Vector4f& viewport,
  const igl::HeatGeodesicsData<double>& data,
  Eigen::VectorXd& D);

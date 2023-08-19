#include <Eigen/Core>
namespace igl {
  class SLIMData;
  class Timer;
  namespace opengl {
    namespace glfw {
      class Viewer;
    }
  }
}

void param_2d_demo_iter(
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F,
    double & uv_scale_param,
    bool & first_iter,
    igl::SLIMData& sData,
    igl::Timer & timer,
    igl::opengl::glfw::Viewer& viewer);


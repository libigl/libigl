#include <igl/read_triangle_mesh.h>
#include <igl/edges.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/avg_edge_length.h>
#include <igl/opengl/glfw/Viewer.h>
#include "tutorial_shared_path.h"
#include <random>

Eigen::MatrixXd V, P;
Eigen::MatrixXi F, E, L;
Eigen::VectorXd R;
double r;

void sample_arrows()
{
  Eigen::RowVector3d minV = V.colwise().minCoeff().array();
  Eigen::RowVector3d maxV = V.colwise().maxCoeff().array();
  std::default_random_engine gen;
  std::uniform_real_distribution<double> dx(minV(0), maxV(0));
  std::uniform_real_distribution<double> dy(minV(1), maxV(1));

  int n = 100;
  P.resize(n * 4, 3);
  R.resize(n * 4);
  L.resize(2 * n, 2);
  for (int i = 0; i < n; ++i) {
    Eigen::RowVector3d h(0, 0, maxV(2) - minV(2));
    Eigen::RowVector3d x(dx(gen), dy(gen), maxV(2) + 0.1 * h(2));
    P.row(4*i+0) = x         ; R(4*i+0) = 0;
    P.row(4*i+1) = x + 0.08*h; R(4*i+1) = 0.3 * r;
    P.row(4*i+2) = x + 0.08*h; R(4*i+2) = 0.1 * r;
    P.row(4*i+3) = x + 0.3*h ; R(4*i+3) = 0.1 * r;
    L.row(2*i+0) << 4*i+0, 4*i+1;
    L.row(2*i+1) << 4*i+2, 4*i+3;
  }
}

int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  // std::string filename = "/home/jdumas/resources/models/dragon.obj";
  std::string filename = TUTORIAL_SHARED_PATH "/cube.obj";
  igl::read_triangle_mesh(filename, V, F);
  igl::edges(F, E);

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.core.align_camera_center(V, F);
  viewer.core.radius_in_screen_space = false;
  Eigen::VectorXd radius(V.rows());
  // int a = E.rows() / 4 + E.rows() / 8;
  // int m = E.rows() / 16;
  // E.topRows(m) = E.middleRows(a, m);
  // E.conservativeResize(m, E.cols());
  // std::cout << E << std::endl;
  r = igl::avg_edge_length(V, F) / 10.0 * viewer.core.camera_base_zoom;
  std::default_random_engine gen;
  std::uniform_real_distribution<double> dist(-1, 1);
  radius = radius.unaryExpr([&](double x) { return r + 0.8 * r * dist(gen); });
  std::cout << radius << std::endl;
  // radius.setConstant(r);
  viewer.data().set_points(V, Eigen::RowVector3d(1, 0, 0), radius);
  viewer.data().set_edges(V, E, Eigen::RowVector3d(1, 0, 0), radius);

  sample_arrows();
  viewer.data().set_edges(P, L, Eigen::RowVector3d(0, 1, 0), R);

  viewer.data().show_overlay = true;
  viewer.data().show_overlay_depth = true;
  viewer.launch();
}

#include <igl/read_triangle_mesh.h>
#include <igl/edges.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/avg_edge_length.h>
#include <igl/opengl/glfw/Viewer.h>
#include "tutorial_shared_path.h"
#include <random>

Eigen::MatrixXd V;
Eigen::MatrixXi F, E;

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
  double r = igl::avg_edge_length(V, F) / 10.0 * viewer.core.camera_base_zoom;
  std::default_random_engine gen;
  std::uniform_real_distribution<double> dist(-1, 1);
  radius = radius.unaryExpr([&](double x) { return r + 0.8 * r * dist(gen); });
  std::cout << radius << std::endl;
  // radius.setConstant(r);
  viewer.data().set_points(V, Eigen::RowVector3d(1, 0, 0), radius);
  viewer.data().set_edges(V, E, Eigen::RowVector3d(1, 0, 0), radius);
  viewer.data().show_overlay = true;
  viewer.data().show_overlay_depth = true;
  viewer.launch();
}

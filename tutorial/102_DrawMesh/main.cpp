#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/bounding_box_diagonal.h>
#include "tutorial_shared_path.h"

Eigen::MatrixXd V;
Eigen::MatrixXi F;

int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  igl::readOFF(TUTORIAL_SHARED_PATH "/bunny.off", V, F);

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  Eigen::VectorXd radius(V.rows());
  radius.setConstant(0.005 / igl::bounding_box_diagonal(V));
  viewer.data().set_points(V, Eigen::RowVector3d(1, 0, 0), radius);
  viewer.launch();
}


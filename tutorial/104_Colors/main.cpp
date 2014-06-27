#define IGL_HEADER_ONLY
#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/jet.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd C;

int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  igl::readOFF("../shared/screwdriver.off", V, F);

  // Plot the mesh
  igl::Viewer viewer;
  viewer.set_mesh(V, F);

  // Use the x coordinate as a scalar field over the surface
  Eigen::VectorXd x = V.col(2);

  // Compute per-vertex colors
  igl::jet(x,true,C);

  // Add per-vertex colors
  viewer.set_colors(C);

  // Launch the viewer
  viewer.launch();
}

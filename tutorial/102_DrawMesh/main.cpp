#define IGL_HEADER_ONLY
#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F;

int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  igl::readOFF("bunny.off", V, F);

  // Plot the mesh
  igl::Viewer viewer;
  viewer.draw_mesh(V, F);
  viewer.launch();
}

#define IGL_HEADER_ONLY
#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/jet.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F;

int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  igl::readOFF("../shared/screwdriver.off", V, F);

  // Plot the mesh
  igl::Viewer viewer;
  viewer.set_mesh(V, F);


  // Normalize x coordinate between 0 and 1
  Eigen::VectorXd value = V.col(0).array() - V.col(0).minCoeff();
  value = value.array() / value.maxCoeff();

  // Map to colors using jet colorramp
  Eigen::MatrixXd C(V.rows(),3);
  for (unsigned i=0; i<V.rows(); ++i)
  {
    double r,g,b;
    igl::jet(value(i),r,g,b);
    C.row(i) << r,g,b;
  }

  // Add per-vertex colors
  viewer.set_colors(C);

  // Launch the viewer
  viewer.launch();
}

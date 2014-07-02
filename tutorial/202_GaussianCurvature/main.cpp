#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/gaussian_curvature.h>
#include <igl/jet.h>

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  MatrixXd V;
  MatrixXi F;
  igl::readOFF("../shared/bumpy.off",V,F);

  VectorXd K;
  igl::gaussian_curvature(V,F,K);

  // Compute pseudocolor
  MatrixXd C;
  igl::jet(K,true,C);

  // Plot the mesh with pseudocolors
  igl::Viewer viewer;
  viewer.data.set_mesh(V, F);
  viewer.data.set_colors(C);
  viewer.launch();
}

#include <igl/viewer/Viewer.h>
#include <igl/tetgen/tetrahedralize.h>
#include <igl/readOFF.h>

// Input polygon
Eigen::MatrixXd V;
Eigen::MatrixXi F;

// Tetrahedralized interior
Eigen::MatrixXd TV;
Eigen::MatrixXi TT;
Eigen::MatrixXi TF;

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;

  // Load a surface mesh
  igl::readOFF("../shared/fertility.off",V,F);
  
  // Tetrahedralize the interior
  igl::tetrahedralize(V,F,"pq1.414", TV,TT,TF);

  // Compute berycenters
  
  // Plot the generated mesh
  igl::Viewer viewer;
  viewer.set_mesh(V,F);
  viewer.launch();
}

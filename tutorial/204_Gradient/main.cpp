#include <igl/readOFF.h>
#include <igl/readDMAT.h>
#include <igl/grad.h>
#include <igl/avg_edge_length.h>
#include <igl/jet.h>
#include <igl/barycenter.h>
#include <igl/viewer/Viewer.h>

#include <iostream>

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  MatrixXd V;
  MatrixXi F;

  // Load a mesh in OFF format
  igl::readOFF("../shared/cheburashka.off", V, F);

  // Read scalar function values from a file, U: #V by 1
  VectorXd U;
  igl::readDMAT("../shared/cheburashka-scalar.dmat",U);

  // Compute gradient operator: #F*3 by #V
  SparseMatrix<double> G;
  igl::grad(V,F,G);

  // Compute gradient of U
  MatrixXd GU = Map<const MatrixXd>((G*U).eval().data(),F.rows(),3);
  // Compute gradient magnitude
  const VectorXd GU_mag = GU.rowwise().norm();

  igl::Viewer viewer;
  viewer.set_mesh(V, F);

  // Compute pseudocolor for original function
  MatrixXd C;
  igl::jet(U,true,C);
  // // Or for gradient magnitude
  //igl::jet(GU_mag,true,C);
  viewer.set_colors(C);

  // Average edge length divided by average gradient (for scaling)
  const double max_size = igl::avg_edge_length(V,F) / GU_mag.mean();
  // Draw a black segment in direction of gradient at face barycenters
  MatrixXd BC;
  igl::barycenter(V,F,BC);
  const RowVector3d black(0,0,0);
  viewer.add_edges(BC,BC+max_size*GU, black);

  // Hide wireframe
  viewer.core.show_lines = false;

  viewer.launch();
}

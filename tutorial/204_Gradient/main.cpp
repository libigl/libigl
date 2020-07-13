#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/grad.h>
#include <igl/jet.h>
#include <igl/readDMAT.h>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>

#include <iostream>
#include "tutorial_shared_path.h"

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  MatrixXd V;
  MatrixXi F;

  // Load a mesh in OFF format
  igl::readOFF(TUTORIAL_SHARED_PATH "/cheburashka.off", V, F);

  // Read scalar function values from a file, U: #V by 1
  VectorXd U;
  igl::readDMAT(TUTORIAL_SHARED_PATH "/cheburashka-scalar.dmat",U);

  // Compute gradient operator: #F*3 by #V
  SparseMatrix<double> G;
  igl::grad(V,F,G);

  // Compute gradient of U
  MatrixXd GU = Map<const MatrixXd>((G*U).eval().data(),F.rows(),3);
  // Compute gradient magnitude
  const VectorXd GU_mag = GU.rowwise().norm();

  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);

  viewer.data().set_data(U);

  // Average edge length divided by average gradient (for scaling)
  const double max_size = igl::avg_edge_length(V,F) / GU_mag.mean();
  // Draw a black segment in direction of gradient at face barycenters
  MatrixXd BC;
  igl::barycenter(V,F,BC);
  const RowVector3d black(0,0,0);
  viewer.data().add_edges(BC,BC+max_size*GU, black);

  // Hide wireframe
  viewer.data().show_lines = false;

  viewer.launch();
}

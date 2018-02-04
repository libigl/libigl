#include <igl/gaussian_curvature.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/jet.h>
#include "tutorial_shared_path.h"

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  MatrixXd V;
  MatrixXi F;
  igl::readOFF(TUTORIAL_SHARED_PATH "/bumpy.off",V,F);

  VectorXd K;
  // Compute integral of Gaussian curvature
  igl::gaussian_curvature(V,F,K);
  // Compute mass matrix
  SparseMatrix<double> M,Minv;
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,M);
  igl::invert_diag(M,Minv);
  // Divide by area to get integral average
  K = (Minv*K).eval();

  // Compute pseudocolor
  MatrixXd C;
  igl::jet(K,true,C);

  // Plot the mesh with pseudocolors
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_colors(C);
  viewer.launch();
}

#include <igl/floor.h>
#include <igl/readOFF.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include "tutorial_shared_path.h"

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  MatrixXd V;
  MatrixXi F;
  igl::readOFF(TUTORIAL_SHARED_PATH "/decimated-knight.off",V,F);

  // 100 random indices into rows of F
  VectorXi I;
  igl::floor((0.5*(VectorXd::Random(100,1).array()+1.)*F.rows()).eval(),I);
  
  // 50 random indices into rows of I
  VectorXi J;
  igl::floor((0.5*(VectorXd::Random(50,1).array()+1.)*I.rows()).eval(),J);
  
  // K = I(J);
  VectorXi K;
  igl::slice(I,J,K);

  // default green for all faces
  MatrixXd C = RowVector3d(0.4,0.8,0.3).replicate(F.rows(),1);
  // Red for each in K
  MatrixXd R = RowVector3d(1.0,0.3,0.3).replicate(K.rows(),1);
  // C(K,:) = R
  igl::slice_into(R,K,1,C);

  // Plot the mesh with pseudocolors
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_colors(C);
  viewer.launch();
}

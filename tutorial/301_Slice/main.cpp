#include <igl/floor.h>
#include <igl/readOFF.h>
#include <igl/find.h>
#include <igl/opengl/glfw/Viewer.h>
#include <iostream>

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

  VectorXi K = I(J);
  // igl::slice(I,J,K); no longer needed

  // default green for all faces
  MatrixXd C = RowVector3d(0.4,0.8,0.3).replicate(F.rows(),1);
  // Red for each in K
  MatrixXd R = RowVector3d(1.0,0.3,0.3).replicate(K.rows(),1);
  // C(K,:) = R
  C(K,Eigen::all) = R;
  // igl::slice_into(R,K,1,C); no longer needed

  Eigen::Array<bool,Eigen::Dynamic,1> W = Eigen::VectorXd::Random(F.rows()).array()>0.5;
  // Set 1/4 of the colors  to blue
  MatrixXd B = RowVector3d(0.3,0.3,1.0).replicate(W.count(),1);
  C(igl::find(W),Eigen::all) = B;

  // Plot the mesh with pseudocolors
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_colors(C);
  viewer.launch();
}

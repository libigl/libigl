
#include <igl/readOFF.h>
#include <igl/cotmatrix.h>
#include <igl/matlab/matlabinterface.h>
#include <igl/opengl/glfw/Viewer.h>
#include <iostream>

// On mac you may need to issue something like:
//
//     PATH=$PATH:/Applications/MATLAB_R2019a.app/bin/ ./tutorial/602_Matlab_bin

// Base mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;

// Matlab instance
Engine* engine;

// Eigenvectors of the laplacian
Eigen::MatrixXd EV;

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  if (key >= '1' && key <= '9')
  {
    viewer.data().set_data(EV.col((key - '1') + 1));
  }

  return false;
}

int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  igl::readOFF(TUTORIAL_SHARED_PATH "/3holes.off", V, F);

  // Launch MATLAB
  igl::matlab::mlinit(&engine);

  // Compute the discrete Laplacian operator
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V,F,L);

  // Send Laplacian matrix to matlab
  igl::matlab::mlsetmatrix(&engine,"L",L);

  // Plot the laplacian matrix using matlab spy
  igl::matlab::mleval(&engine,"spy(L)");

  // Extract the first 10 eigenvectors
  igl::matlab::mleval(&engine,"[EV,~] = eigs(-L,10,'sm')");

  // Plot the size of EV (only for demonstration purposes)
  std::cerr << igl::matlab::mleval(&engine,"size(EV)") << std::endl;

  // Retrieve the result
  igl::matlab::mlgetmatrix(&engine,"EV",EV);

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.callback_key_down = &key_down;
  viewer.data().set_mesh(V, F);

  // Plot the first non-trivial eigenvector
  viewer.data().set_data(EV.col(1));

  // Launch the viewer
  viewer.launch();
}

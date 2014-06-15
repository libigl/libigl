#define IGL_HEADER_ONLY
#include <igl/readOFF.h>
#include <igl/writeOBJ.h>
#include <iostream>

Eigen::MatrixXd V;
Eigen::MatrixXi F;

int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  igl::readOFF("cube.off", V, F);

  // Plot the vertices and faces matrices
  std::cerr << "Vertices: " << std::endl << V << std::endl;
  std::cerr << "Faces:    " << std::endl << F << std::endl;

  // Save the mesh in OBJ format
  igl::writeOBJ("cube.obj",V,F);
}

#include <igl/readOFF.h>
#include <iostream>

Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> V;
Eigen::MatrixXi F;

int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  igl::readOFF(TUTORIAL_SHARED_PATH  "/cube.off", V, F);
}

#include <Eigen/Core>

#include <igl/readDMAT.h>
#include <cstdio>

#ifdef IGL_STATIC_LIBRARY
#  define IGL_HEADER_ONLY
#  define IGL_HEADER_ONLY_WAS_NOT_DEFINED
#endif
#include <igl/matlab/MatlabWorkspace.h>
#include <igl/on_boundary.h>
#ifndef IGL_STATIC_LIBRARY_WAS_NOT_DEFINED
#  undef IGL_HEADER_ONLY
#endif

int main(int argc, char * argv[])
{
  using namespace igl;
  using namespace Eigen;
  if(argc <= 2)
  {
    printf("USAGE:\n  ./example [input.dmat] [output.mat]\n");
    return 1;
  }
  MatrixXd M;
  readDMAT(argv[1],M);
  igl::matlab::MatlabWorkspace mat;
  mat.save(M,"M");
  mat.write(argv[2]);

  return 0;
}

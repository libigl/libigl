#include <Eigen/Core>

#include <igl/readDMAT.h>
#include <cstdio>

#ifndef IGL_HEADER_ONLY
#  define IGL_HEADER_ONLY
#  define IGL_HEADER_ONLY_WAS_NOT_DEFINED
#endif
#include <igl/matlab/MatlabWorkspace.h>
#include <igl/on_boundary.h>
#ifdef IGL_HEADER_ONLY_WAS_NOT_DEFINED
#  undef IGL_HEADER_ONLY
#endif

int main(int argc, char * argv[])
{
  using namespace igl;
  using namespace Eigen;
  if(argc <= 2)
  {
    printf("USAGE:\n  ./example [input path] [output path]\n");
    return 1;
  }
  MatrixXd M;
  readDMAT(argv[1],M);
  MatlabWorkspace mat;
  mat.save(M,"M");
  mat.write(argv[2]);

  return 0;
}

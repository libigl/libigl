#include <Eigen/Core>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
using namespace igl;
#include <cstdio>

int main(int argc, char * argv[])
{
  if(argc <= 2)
  {
    printf("USAGE:\n  ./example [input path] [output path]\n");
    return 1;
  }
  Eigen::MatrixXd M;
  readDMAT(argv[1],M);
  writeDMAT(argv[2],M);
  return 0;
}

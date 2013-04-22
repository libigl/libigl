// Compile with GCC:
// g++ -I/usr/local/igl/igl_lib/include -DIGL_HEADER_ONLY \
//   -I/opt/local/include/eigen3 example.cpp -o example
// Run with:
// ./example ../shared/TinyTorus.obj
#include <Eigen/Dense>
#include <igl/readOBJ.h>
#include <iostream>
int main(int argc, char * argv[])
{
  if(argc>1)
  {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::readOBJ(argv[1],V,F);
    std::cout<<"Hello, mesh with "<<V.rows()<<" vertices!"<<std::endl;
  }else{
    std::cout<<"Hello, world!"<<std::endl;
  }
  return 0;
}


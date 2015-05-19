#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <string>
#include <iostream>
int main(int argc, char * argv[])
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  MatrixXd V;
  MatrixXi F;
  string in,out;
  switch(argc)
  {
    case 3:
      in = argv[1];
      out = argv[2];
      break;
    default:
      cerr<<R"(
USAGE:
  convertmesh input.[mesh|obj|off|ply|stl|wrl] output.[mesh|obj|off|ply|stl|wrl]

  Note: .ply and .stl outputs are binary.
)";
    return EXIT_FAILURE;
  }
  return 
    read_triangle_mesh(in,V,F) && write_triangle_mesh(out,V,F,false) ? 
    EXIT_SUCCESS : EXIT_FAILURE;
}

#include <igl/opengl/glfw/Viewer.h>
#include <igl/triangle/triangulate.h>
// Input polygon
Eigen::MatrixXd V;
Eigen::MatrixXi E;
Eigen::MatrixXd H;

// Triangulated interior
Eigen::MatrixXd V2;
Eigen::MatrixXi F2;

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;

  // Create the boundary of a square
  V.resize(8,2);
  E.resize(8,2);
  H.resize(1,2);

  V << -1,-1, 1,-1, 1,1, -1, 1,
       -2,-2, 2,-2, 2,2, -2, 2;

  E << 0,1, 1,2, 2,3, 3,0,
       4,5, 5,6, 6,7, 7,4;

  H << 0,0;

  // Triangulate the interior
  igl::triangle::triangulate(V,E,H,"a0.005q",V2,F2);

  // Plot the generated mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V2,F2);
  viewer.launch();
}

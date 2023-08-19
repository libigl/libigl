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

  // create two squares, one with edge length of 4,
  // one with edge length of 2
  // both centered at origin
  V << -1,-1, 1,-1, 1,1, -1, 1,
       -2,-2, 2,-2, 2,2, -2, 2;

  // add the edges of the squares
  E << 0,1, 1,2, 2,3, 3,0,
       4,5, 5,6, 6,7, 7,4;

  // specify a point that is inside a closed shape
  // where we do not want triangulation to happen
  H << 0,0;

  // Triangulate the interior
  // a0.005 means that the area of each triangle should
  // not be greater than 0.005
  // q means that no angles will be smaller than 20 degrees
  // for a detailed set of commands please refer to:
  // https://www.cs.cmu.edu/~quake/triangle.switch.html
  igl::triangle::triangulate(V,E,H,"a0.005q",V2,F2);

  // Plot the generated mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V2,F2);
  viewer.launch();
}

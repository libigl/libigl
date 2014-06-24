#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/local_basis.h>
#include <igl/barycenter.h>
#include <igl/avg_edge_length.h>
#include <igl/embree/ambient_occlusion.h>
#include <algorithm>

using namespace Eigen;
using namespace std;

// Mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;

VectorXd AO;

  // It allows to change the degree of the field when a number is pressed
bool key_down(igl::Viewer& viewer, unsigned char key, int modifier)
{
  if (key == '1')
  {
    // Show the mesh without the ambient occlusion factor
    viewer.set_colors(Eigen::RowVector3d(1,1,1));
  }

  if (key == '2')
  {
    // Show the mesh with the ambient occlusion factor
    MatrixXd C = MatrixXd::Constant(V.rows(),3,1);
    for (unsigned i=0; i<C.rows();++i)
      C.row(i) *= std::min<double>(AO(i)+0.2,1);
    viewer.set_colors(C);
  }

  return false;
}

int main(int argc, char *argv[])
{
  using namespace std;
  using namespace Eigen;

  // Load a mesh in OFF format
  igl::readOFF("../shared/fertility.off", V, F);

  MatrixXd N;
  igl::per_vertex_normals(V,F,N);

  // Compute ambient occlusion factor using embree
  igl::ambient_occlusion(V,F,V,N,500,AO);
  AO = 1.0 - AO.array();

  // Show mesh
  igl::Viewer viewer;
  viewer.set_mesh(V, F);
  viewer.callback_key_down = &key_down;
  key_down(viewer,'2',0);
  viewer.options.show_lines = false;
  viewer.launch();
}

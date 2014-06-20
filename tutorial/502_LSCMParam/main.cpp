#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/boundary_vertices_sorted.h>

#include <igl/lscm.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd V_uv;

bool key_down(igl::Viewer& viewer, unsigned char key, int modifier)
{

  if (key == '1')
    // Plot the 3D mesh
    viewer.set_mesh(V,F);
  else if (key == '2')
    // Plot the mesh in 2D using the UV coordinates as vertex coordinates
    viewer.set_mesh(V_uv,F);

  viewer.compute_normals();

  return false;
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;

  // Load a mesh in OFF format
  igl::readOFF("../shared/camelhead.off", V, F);

  // Fix two points on the boundary
  VectorXi bnd,b(2,1);
  igl::boundary_vertices_sorted(V,F,bnd);
  b(0) = bnd(0);
  b(1) = bnd(round(bnd.size()/2));
  MatrixXd bc(2,2);
  bc<<0,0,1,0;

  // LSCM parametrization
  igl::lscm(V,F,b,bc,V_uv);

  // Scale the uv
  V_uv *= 5;

  // Plot the mesh
  igl::Viewer viewer;
  viewer.set_mesh(V, F);
  viewer.set_uv(V_uv);
  viewer.callback_key_down = &key_down;

  // Disable wireframe
  viewer.options.show_lines = false;

  // Draw checkerboard texture
  viewer.options.show_texture = true;

  // Launch the viewer
  viewer.launch();
}

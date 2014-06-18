#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/boundary_vertices_sorted.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>

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
  // Load a mesh in OFF format
  igl::readOFF("../shared/camelhead.off", V, F);

  // Find the open boundary
  Eigen::VectorXi bnd;
  igl::boundary_vertices_sorted(V,F,bnd);

  // Map the boundary to a circle, preserving edge proportions
  Eigen::MatrixXd bnd_uv;
  igl::map_vertices_to_circle(V,F,bnd,bnd_uv);

  // Harmonic parametrization for the internal vertices
  igl::harmonic(V,F,bnd,bnd_uv,1,V_uv);

  // Scale UV to make the texture more clear
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

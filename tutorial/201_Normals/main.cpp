#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_corner_normals.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F;

Eigen::MatrixXd N_vertices;
Eigen::MatrixXd N_faces;
Eigen::MatrixXd N_corners;


// This function is called every time a keyboard button is pressed
bool key_down(igl::Viewer& viewer, unsigned char key, int modifier)
{
  switch(key)
  {
    case '1':
      viewer.set_normals(N_faces);
      return true;
    case '2':
      viewer.set_normals(N_vertices);
      return true;
    case '3':
      viewer.set_normals(N_corners);
      return true;
    default: break;
  }
  return false;
}

int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  igl::readOFF("../shared/fandisk.off", V, F);

  // Compute per-face normals
  igl::per_face_normals(V,F,N_faces);

  // Compute per-vertex normals
  igl::per_vertex_normals(V,F,N_vertices);

  // Compute per-corner normals, |dihedral angle| > 20 degrees --> crease
  igl::per_corner_normals(V,F,20,N_corners);

  // Plot the mesh
  igl::Viewer viewer;
  viewer.callback_key_down = &key_down;
  viewer.core.show_lines = false;
  viewer.set_mesh(V, F);
  viewer.set_normals(N_vertices);
  viewer.launch();
}

#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_corner_normals.h>
#include <iostream>
#include "tutorial_shared_path.h"

Eigen::MatrixXd V;
Eigen::MatrixXi F;

Eigen::MatrixXd N_vertices;
Eigen::MatrixXd N_faces;
Eigen::MatrixXd N_corners;


// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  switch(key)
  {
    case '1':
      viewer.data().set_normals(N_faces);
      return true;
    case '2':
      viewer.data().set_normals(N_vertices);
      return true;
    case '3':
      viewer.data().set_normals(N_corners);
      return true;
    default: break;
  }
  return false;
}

int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  igl::read_triangle_mesh(
    argc>1?argv[1]: TUTORIAL_SHARED_PATH "/fandisk.off",V,F);

  // Compute per-face normals
  igl::per_face_normals(V,F,N_faces);

  // Compute per-vertex normals
  igl::per_vertex_normals(V,F,N_vertices);

  // Compute per-corner normals, |dihedral angle| > 20 degrees --> crease
  igl::per_corner_normals(V,F,20,N_corners);

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.callback_key_down = &key_down;
  viewer.data().show_lines = false;
  viewer.data().set_mesh(V, F);
  viewer.data().set_normals(N_faces);
  std::cout<<
    "Press '1' for per-face normals."<<std::endl<<
    "Press '2' for per-vertex normals."<<std::endl<<
    "Press '3' for per-corner normals."<<std::endl;
  viewer.launch();
}

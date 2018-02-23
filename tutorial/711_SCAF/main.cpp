#include <igl/scaf.h>
#include <igl/arap.h>
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>

#include "tutorial_shared_path.h"

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd V_uv;
Eigen::MatrixXd initial_guess;
igl::SCAFData scaf_data;

bool show_uv = false;
float uv_scale = 0.2;

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  if (key == '1')
    show_uv = false;
  else if (key == '2')
    show_uv = true;

  if (key == 'q')
    V_uv = initial_guess;
  if (key == ' ')
    igl::scaf_solve(scaf_data, 1);

  const auto& V_uv = uv_scale * scaf_data.w_uv.topRows(V.rows());
  if (show_uv)
  {
    viewer.data().clear();
    viewer.data().set_mesh(V_uv,F);
    viewer.data().set_uv(V_uv);
    viewer.core.align_camera_center(V_uv,F);
  }
  else
  {
    viewer.data().set_mesh(V,F);
    viewer.data().set_uv(V_uv);
    viewer.core.align_camera_center(V,F);
  }

  viewer.data().compute_normals();

  return false;
}

int main(int argc, char *argv[])
{
  using namespace std;
  // Load a mesh in OFF format
  igl::readOBJ(TUTORIAL_SHARED_PATH "/camel_b.obj", V, F);

  Eigen::VectorXi b;
  Eigen::MatrixXd bc, V_init;
  igl::scaf_precompute(V, F, V_init, scaf_data, igl::SLIMData::SYMMETRIC_DIRICHLET, b, bc, 0);

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  const auto& V_uv = uv_scale * scaf_data.w_uv.topRows(V.rows());
  viewer.data().set_uv(V_uv);
  viewer.callback_key_down = &key_down;

  // Enable wireframe
  viewer.data().show_lines = true;

  // Draw checkerboard texture
  viewer.data().show_texture = true;

  // Launch the viewer
  viewer.launch();
}

#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/svd3x3/arap.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd V_uv;
Eigen::MatrixXd initial_guess;

bool show_uv = false;

bool key_down(igl::Viewer& viewer, unsigned char key, int modifier)
{
  if (key == '1')
    show_uv = false;
  else if (key == '2')
    show_uv = true;

  if (key == 'q')
    V_uv = initial_guess;

  if (show_uv)
    viewer.set_mesh(V_uv,F);
  else
    viewer.set_mesh(V,F);

  viewer.compute_normals();

  return false;
}

int main(int argc, char *argv[])
{
  using namespace std;
  // Load a mesh in OFF format
  igl::readOFF("../shared/camelhead.off", V, F);

  // Compute the initial solution for ARAP (harmonic parametrization)
  Eigen::VectorXi bnd;
  igl::boundary_loop(V,F,bnd);
  Eigen::MatrixXd bnd_uv;
  igl::map_vertices_to_circle(V,bnd,bnd_uv);

  igl::harmonic(V,F,bnd,bnd_uv,1,initial_guess);

  // Add dynamic regularization to avoid to specify boundary conditions
  igl::ARAPData arap_data;
  arap_data.with_dynamics = true;
  Eigen::VectorXi b  = Eigen::VectorXi::Zero(0);
  Eigen::MatrixXd bc = Eigen::MatrixXd::Zero(0,0);
  
  // Initialize ARAP
  arap_data.max_iter = 100;
  // 2 means that we're going to *solve* in 2d
  arap_precomputation(V,F,2,b,arap_data);


  // Solve arap using the harmonic map as initial guess
  V_uv = initial_guess;

  arap_solve(bc,arap_data,V_uv);


  // Scale UV to make the texture more clear
  V_uv *= 20;

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

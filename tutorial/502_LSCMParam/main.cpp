#include <igl/boundary_loop.h>
#include <igl/read_triangle_mesh.h>
#include <igl/viewer/Viewer.h>

#include <igl/lscm.h>

#include "tutorial_shared_path.h"
#include <iostream>


Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd V_uv;

bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{

  if (key == '1')
  {
    // Plot the 3D mesh
    viewer.data.set_mesh(V,F);
    viewer.data.compute_normals();
    viewer.core.align_camera_center(V,F);
  }
  else if (key == '2')
  {
    // Plot the mesh in 2D using the UV coordinates as vertex coordinates
    viewer.data.set_mesh(V_uv,F);
    viewer.data.compute_normals();
    viewer.core.align_camera_center(V_uv,F);
  }


  return false;
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;

  // Set input mesh filename
  std::string filename(TUTORIAL_SHARED_PATH "/camelhead.off");
  if (argc > 1)
      filename = std::string(argv[1]);

  // Try to load the input mesh
  if (igl::read_triangle_mesh(filename, V, F) == false)
      return -1;

  // Fix two points on the boundary
  VectorXi bnd,b(2,1);
  igl::boundary_loop(F,bnd);
  if (bnd.size() < 1)
        std::cerr << "error: Mesh has no boundary"<<std::endl;
  b(0) = bnd(0);
  b(1) = bnd(round(bnd.size()/2));
  MatrixXd bc(2,2);
  bc<<0,0,1,0;

  // LSCM parametrization
  igl::lscm(V,F,b,bc,V_uv);

  // Scale the uv
  V_uv *= 5;

  // Plot the mesh
  igl::viewer::Viewer viewer;
  viewer.data.set_mesh(V, F);
  viewer.data.set_uv(V_uv);
  viewer.callback_key_down = &key_down;

  // Disable wireframe
  viewer.core.show_lines = false;

  // Draw checkerboard texture
  viewer.core.show_texture = true;

  // Launch the viewer
  viewer.launch();
}

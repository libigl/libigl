#include "tutorial_shared_path.h"
#include <igl/readOFF.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <iostream>

int main(int argc, char *argv[])
{
  // Mesh with per-face color
  Eigen::MatrixXd V, C;
  Eigen::MatrixXi F;

  // Load a mesh in OFF format
  igl::readOFF(TUTORIAL_SHARED_PATH "/fertility.off", V, F);

  // Initialize white
  C = Eigen::MatrixXd::Constant(F.rows(),3,1);
  igl::opengl::glfw::Viewer viewer;
  viewer.callback_mouse_down =
    [&V,&F,&C](igl::opengl::glfw::Viewer& viewer, int, int)->bool
  {
    int fid;
    Eigen::Vector3f bc;
    // Cast a ray in the view direction starting from the mouse position
    double x = viewer.current_mouse_x;
    double y = viewer.core.viewport(3) - viewer.current_mouse_y;
    if(igl::unproject_onto_mesh(Eigen::Vector2f(x,y), viewer.core.view,
      viewer.core.proj, viewer.core.viewport, V, F, fid, bc))
    {
      // paint hit red
      C.row(fid)<<1,0,0;
      viewer.data().set_colors(C);
      return true;
    }
    return false;
  };
  std::cout<<R"(Usage:
  [click]  Pick face on shape

)";
  // Show mesh
  viewer.data().set_mesh(V, F);
  viewer.data().set_colors(C);
  viewer.data().show_lines = false;
  viewer.launch();
}

#include "tutorial_shared_path.h"
#include <igl/read_triangle_mesh.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/viewer/Viewer.h>
#include <iostream>

int main(int argc, char *argv[])
{
  // Mesh with per-face color
  Eigen::MatrixXd V, C;
  Eigen::MatrixXi F;

  // Set input mesh filename
  std::string filename(TUTORIAL_SHARED_PATH "/fertility.off");
  if (argc > 1)
      filename = std::string(argv[1]);

  // Try to load the input mesh
  if (igl::read_triangle_mesh(filename, V, F) == false)
      return -1;

  // Initialize white
  C = Eigen::MatrixXd::Constant(F.rows(),3,1);
  igl::viewer::Viewer viewer;
  viewer.callback_mouse_down = 
    [&V,&F,&C](igl::viewer::Viewer& viewer, int, int)->bool
  {
    int fid;
    Eigen::Vector3f bc;
    // Cast a ray in the view direction starting from the mouse position
    double x = viewer.current_mouse_x;
    double y = viewer.core.viewport(3) - viewer.current_mouse_y;
    if(igl::unproject_onto_mesh(Eigen::Vector2f(x,y), viewer.core.view * viewer.core.model,
      viewer.core.proj, viewer.core.viewport, V, F, fid, bc))
    {
      // paint hit red
      C.row(fid)<<1,0,0;
      viewer.data.set_colors(C);
      std::cout << "Picked: " << fid << std::endl;
      return true;
    }
    return false;
  };
  std::cout<<R"(Usage:
  [click]  Pick face on shape

)";
  // Show mesh
  viewer.data.set_mesh(V, F);
  viewer.data.set_colors(C);
  viewer.core.show_lines = false;
  viewer.launch();
}

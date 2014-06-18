#define IGL_HEADER_ONLY
#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>

Eigen::MatrixXd V1,V2;
Eigen::MatrixXi F1,F2;

// This function is called every time a keyboard button is pressed
bool key_down(igl::Viewer& viewer, unsigned char key, int modifier)
{
  if (key == '1')
  {
    // Clear mesh should be called before drawing the mesh
    viewer.clear_mesh();
    // Draw_mesh creates or updates the vertices and faces of the displayed mesh.
    // If a mesh is already displayed, draw_mesh returns an error if the given V and
    // F have size different than the current ones
    viewer.set_mesh(V1, F1);
  }
  else if (key == '2')
  {
    viewer.clear_mesh();
    viewer.set_mesh(V2, F2);
  }

  return false;
}


int main(int argc, char *argv[])
{
  // Load two meshes
  igl::readOFF("../shared/bumpy.off", V1, F1);
  igl::readOFF("../shared/fertility.off", V2, F2);

  igl::Viewer viewer;
  // Register a keyboard callback that allows to switch between
  // the two loaded meshes
  viewer.callback_key_down = &key_down;
  viewer.set_mesh(V1, F1);
  viewer.launch();
}

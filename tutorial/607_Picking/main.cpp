#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/embree/EmbreeIntersector.h>
#include <igl/embree/unproject_onto_mesh.h>

#include <algorithm>

#include "tutorial_shared_path.h"

// Mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;

// Per-vertex color
Eigen::MatrixXd C;

igl::embree::EmbreeIntersector* ei;

std::vector<int> picked_vertices;

bool mouse_down(igl::viewer::Viewer& viewer, int button, int modifier)
{

  using namespace Eigen;
  using namespace std;

  int vid, fid;

  // Cast a ray in the view direction starting from the mouse position
  double x = viewer.current_mouse_x;
  double y = viewer.core.viewport(3) - viewer.current_mouse_y;
  bool hit = igl::embree::unproject_onto_mesh(Vector2f(x,y),
                                F,
                                viewer.core.view * viewer.core.model,
                                viewer.core.proj,
                                viewer.core.viewport,
                                *ei,
                                fid,
                                vid);

  // If the ray hits a face, assign a red color to the closest vertex
  if (hit)
  {
    cerr << "Picked face(vertex): " << fid << " (" << vid << ")" << endl;
    C.row(vid) << 1,0,0;
    viewer.data.set_colors(C);
    return true;
  }

  return false;
}

int main(int argc, char *argv[])
{

  using namespace Eigen;
  using namespace std;
  // Load a mesh in OFF format
  igl::readOFF(TUTORIAL_SHARED_PATH "/fertility.off", V, F);

  // Create a BVH for raycasting
  ei = new igl::embree::EmbreeIntersector();
  ei->init(V.cast<float>(),F);

  // Initialize the colors to white for all vertices
  C = MatrixXd::Constant(V.rows(),3,1);

  // Show mesh
  igl::viewer::Viewer viewer;
  viewer.data.set_mesh(V, F);
  viewer.callback_mouse_down = &mouse_down;
  viewer.data.set_colors(C);
  viewer.core.show_lines = false;
  viewer.launch();
}

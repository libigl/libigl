#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/embree/EmbreeIntersector.h>
#include <igl/embree/unproject_in_mesh.h>

#include <algorithm>

using namespace Eigen;
using namespace std;

// Mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;

// Per-vertex color
MatrixXd C;

igl::EmbreeIntersector* ei;

vector<int> picked_vertices;

bool mouse_down(igl::Viewer& viewer, int button, int modifier)
{

  int vid, fid;

  // Cast a ray in the view direction starting from the mouse position
  double x = viewer.current_mouse_x;
  double y = viewer.viewport(3) - viewer.current_mouse_y;
  bool hit = unproject_in_mesh(Vector2f(x,y),
                                F,
                                viewer.view * viewer.model,
                                viewer.proj,
                                viewer.viewport,
                                *ei,
                                fid,
                                vid);

  // If the ray hits a face, assign a red color to the closest vertex
  if (hit)
  {
    cerr << "Picked face(vertex): " << fid << " (" << vid << ")" << endl;
    C.row(vid) << 1,0,0;
    viewer.set_colors(C);
    return true;
  }

  return false;
}

int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  igl::readOFF("../shared/fertility.off", V, F);

  // Create a BVH for raycasting
  ei = new igl::EmbreeIntersector();
  ei->init(V.cast<float>(),F);

  // Initialize the colors to white for all vertices
  C = MatrixXd::Constant(V.rows(),3,1);

  // Show mesh
  igl::Viewer viewer;
  viewer.set_mesh(V, F);
  viewer.callback_mouse_down = &mouse_down;
  viewer.set_colors(C);
  viewer.core.show_lines = false;
  viewer.launch();
}

#define IGL_HEADER_ONLY
#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <sstream>

Eigen::MatrixXd V;
Eigen::MatrixXi F;

int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  igl::readOFF("../shared/bunny.off", V, F);

  // Find the bounding box
  Eigen::Vector3d m = V.colwise().minCoeff();
  Eigen::Vector3d M = V.colwise().maxCoeff();

  // Corners of the bounding box
  Eigen::MatrixXd V_box(8,3);
  V_box <<
  m(0), m(1), m(2),
  M(0), m(1), m(2),
  M(0), M(1), m(2),
  m(0), M(1), m(2),
  m(0), m(1), M(2),
  M(0), m(1), M(2),
  M(0), M(1), M(2),
  m(0), M(1), M(2);

  // Edges of the bounding box
  Eigen::MatrixXi E_box(12,2);
  E_box <<
  0, 1,
  1, 2,
  2, 3,
  3, 0,
  4, 5,
  5, 6,
  6, 7,
  7, 4,
  0, 4,
  1, 5,
  2, 6,
  7 ,3;

  // Plot the mesh
  igl::Viewer viewer;
  viewer.set_mesh(V, F);

  // Plot the corners of the bounding box as points
  viewer.add_points(V_box,Eigen::RowVector3d(1,0,0));

  // Plot the edges of the bounding box
  for (unsigned i=0;i<E_box.rows(); ++i)
    viewer.add_edges
    (
      V_box.row(E_box(i,0)),
      V_box.row(E_box(i,1)),
      Eigen::RowVector3d(1,0,0)
    );

  // Increase the thickness of the lines
  viewer.options.line_width = 2.0f;

  // Plot labels with the coordinates of bounding box vertices
  std::stringstream l1;
  l1 << m(0) << ", " << m(1) << ", " << m(2);
  viewer.add_label(m,l1.str());
  std::stringstream l2;
  l2 << M(0) << ", " << M(1) << ", " << M(2);
  viewer.add_label(M,l2.str());

  // Launch the viewer
  viewer.launch();
}

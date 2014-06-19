#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_corner_normals.h>
#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/principal_curvature.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F;

int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  igl::readOFF("../shared/fertility.off", V, F);

  // Compute curvature directions via quadric fitting
  Eigen::MatrixXd PD1,PD2,PV1,PV2;
  igl::principal_curvature(V,F,PD1,PD2,PV1,PV2);

  // Plot the mesh
  igl::Viewer viewer;
  viewer.set_mesh(V, F);

  // Find the average edge length
  double avg = igl::avg_edge_length(V,F);
  
  // Draw a red segment on each vertex parallel to the minimal curvature direction
  viewer.add_edges(V + PD1*avg, V - PD1*avg, Eigen::RowVector3d(0,0,1));

  // Draw a blue segment on each vertex parallel to the maximal curvature direction
  viewer.add_edges(V + PD2*avg, V - PD2*avg, Eigen::RowVector3d(1,0,0));

  // Launch the viewer
  viewer.launch();
}

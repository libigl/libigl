#include <igl/opengl/glfw/Viewer.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/readOBJ.h>
#include <igl/marching_tets.h>
#include <Eigen/Core>

#include "tutorial_shared_path.h"


int main(int argc, char * argv[])
{

  // Load a surface mesh which is a cube
  Eigen::MatrixXd surfaceV;
  Eigen::MatrixXi surfaceF;
  igl::readOBJ(TUTORIAL_SHARED_PATH "/cube.obj", surfaceV, surfaceF);

  // Find the centroid of the loaded mesh
  Eigen::RowVector3d surfaceCenter = surfaceV.colwise().sum() / surfaceV.rows();

  // Center the mesh about the origin
  surfaceV.rowwise() -= surfaceCenter;

  // Tetrahedralize the surface mesh
  Eigen::MatrixXd TV; // Tet mesh vertices
  Eigen::MatrixXi TF; // Tet mesh boundary face indices
  Eigen::MatrixXi TT; // Tet mesh tetrahedron indices
  igl::copyleft::tetgen::tetrahedralize(surfaceV, surfaceF, "pq1.414a0.0001", TV, TT, TF);

  // Compute a scalar at each tet vertex which is the distance from the vertex to the origin
  Eigen::VectorXd S = TV.rowwise().norm();

  // Compute a mesh (stored in SV, SF) representing the iso-level-set for the isovalue 0.5
  Eigen::MatrixXd SV;
  Eigen::MatrixXi SF;
  igl::marching_tets(TV, TT, S, 0.45, SV, SF);

  // Draw the mesh stored in (SV, SF)
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(SV, SF);
  viewer.callback_key_down =
    [&](igl::opengl::glfw::Viewer & viewer, unsigned char key, int mod)->bool
    {
      viewer.data().set_face_based(true);
      return true;
    };
  viewer.launch();
}

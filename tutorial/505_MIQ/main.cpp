#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/comiso/miq.h>
#include <igl/barycenter.h>
#include <igl/avg_edge_length.h>
#include <igl/comiso/nrosy.h>
#include <sstream>
#include <igl/rotate_vectors.h>


Eigen::VectorXi Seams;

// Cuts
Eigen::VectorXi C;

// Singularities
Eigen::VectorXd S;

// Cross field
Eigen::MatrixXd X;
Eigen::MatrixXd X2;

// Create a texture that hides the integer translation in the parametrization
void line_texture(Eigen::Matrix<char,Eigen::Dynamic,Eigen::Dynamic> &texture_R,
                  Eigen::Matrix<char,Eigen::Dynamic,Eigen::Dynamic> &texture_G,
                  Eigen::Matrix<char,Eigen::Dynamic,Eigen::Dynamic> &texture_B)
  {
    unsigned size = 128;
    unsigned size2 = size/2;
    unsigned lineWidth = 3;
    texture_R.setConstant(size, size, 255);
    for (unsigned i=0; i<size; ++i)
      for (unsigned j=size2-lineWidth; j<=size2+lineWidth; ++j)
        texture_R(i,j) = 0;
    for (unsigned i=size2-lineWidth; i<=size2+lineWidth; ++i)
      for (unsigned j=0; j<size; ++j)
        texture_R(i,j) = 0;

    texture_G = texture_R;
    texture_B = texture_R;
  }

int main(int argc, char *argv[])
{
  using namespace Eigen;

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  // Load a mesh in OFF format
  igl::readOFF("../shared/3holes.off", V, F);

  // Contrain one face
  VectorXi b(1);
  b << 0;
  MatrixXd bc(1,3);
  bc << 1, 0, 0;

  // Create a smooth 4-RoSy field
  igl::nrosy(V,F,b,bc,VectorXi(),VectorXd(),MatrixXd(),4,0.5,X,S);

  // Find the the orthogonal vector
  MatrixXd B1,B2,B3;
  igl::local_basis(V,F,B1,B2,B3);
  X2 = igl::rotate_vectors(X, VectorXd::Constant(1,M_PI/2), B1, B2);
  
  Eigen::MatrixXd UV;
  Eigen::MatrixXi FUV;

  double gradient_size = 50;
  double iter = 0;
  double stiffness = 5.0;
  bool direct_round = 0;
  igl::miq(V,
           F,
           X,
           X2,
           UV,
           FUV,
           gradient_size,
           stiffness,
           direct_round,
           iter);


  // Face barycenters
  Eigen::MatrixXd MF;
  igl::barycenter(V, F, MF);

  double scale =  .5*igl::avg_edge_length(V, F);

  // Plot the mesh
  igl::Viewer viewer;
  viewer.set_mesh(V, F);

  // Plot the field
  viewer.add_edges (MF, MF+scale*X ,Eigen::RowVector3d(1,0,1));
  viewer.add_edges (MF, MF+scale*X2,Eigen::RowVector3d(1,0,1));
  viewer.set_uv(UV,FUV);
  viewer.options.show_texture = true;

  Eigen::Matrix<char,Eigen::Dynamic,Eigen::Dynamic> texture_R, texture_G, texture_B;
  line_texture(texture_R, texture_G, texture_B);
  viewer.set_texture(texture_R, texture_B, texture_G);
  // Increase the thickness of the lines
  viewer.options.line_width = 2.0f;

  // Launch the viewer
  viewer.launch();
}

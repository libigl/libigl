#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/comb_cross_field.h>
#include <igl/comb_frame_field.h>
#include <igl/compute_frame_field_bisectors.h>
#include <igl/cross_field_missmatch.h>
#include <igl/cut_mesh_from_singularities.h>
#include <igl/find_cross_field_singularities.h>
#include <igl/local_basis.h>
#include <igl/readOFF.h>
#include <igl/rotate_vectors.h>
#include <igl/copyleft/comiso/miq.h>
#include <igl/copyleft/comiso/nrosy.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/PI.h>
#include <sstream>

#include "tutorial_shared_path.h"
#include <igl/serialize.h>

// Input mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;


// Face barycenters
Eigen::MatrixXd B;

// Scale for visualizing the fields
double global_scale;
bool extend_arrows = false;

// Cross field
Eigen::MatrixXd X1,X2;

// Bisector field
Eigen::MatrixXd BIS1, BIS2;

// Combed bisector
Eigen::MatrixXd BIS1_combed, BIS2_combed;

// Per-corner, integer mismatches
Eigen::Matrix<int, Eigen::Dynamic, 3> MMatch;

// Field singularities
Eigen::Matrix<int, Eigen::Dynamic, 1> isSingularity, singularityIndex;

// Per corner seams
Eigen::Matrix<int, Eigen::Dynamic, 3> Seams;

// Combed field
Eigen::MatrixXd X1_combed, X2_combed;


// Global parametrization (with seams)
Eigen::MatrixXd UV_seams;
Eigen::MatrixXi FUV_seams;

// Global parametrization
Eigen::MatrixXd UV;
Eigen::MatrixXi FUV;


// Create a texture that hides the integer translation in the parametrization
void line_texture(Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> &texture_R,
                  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> &texture_G,
                  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> &texture_B)
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

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  if (key == 'E')
  {
    extend_arrows = !extend_arrows;
  }

  if (key <'1' || key >'8')
    return false;

  viewer.data().clear();
  viewer.data().show_lines = false;
  viewer.data().show_texture = false;

  if (key == '1')
  {
    // Cross field
    viewer.data().set_mesh(V, F);
    viewer.data().add_edges(extend_arrows ? B - global_scale*X1 : B, B + global_scale*X1 ,Eigen::RowVector3d(1,0,0));
    viewer.data().add_edges(extend_arrows ? B - global_scale*X2 : B, B + global_scale*X2 ,Eigen::RowVector3d(0,0,1));
  }

  if (key == '2')
  {
    // Bisector field
    viewer.data().set_mesh(V, F);
    viewer.data().add_edges(extend_arrows ? B - global_scale*BIS1 : B, B + global_scale*BIS1 ,Eigen::RowVector3d(1,0,0));
    viewer.data().add_edges(extend_arrows ? B - global_scale*BIS2 : B, B + global_scale*BIS2 ,Eigen::RowVector3d(0,0,1));
  }

  if (key == '3')
  {
    // Bisector field combed
    viewer.data().set_mesh(V, F);
    viewer.data().add_edges(extend_arrows ? B - global_scale*BIS1_combed : B, B + global_scale*BIS1_combed ,Eigen::RowVector3d(1,0,0));
    viewer.data().add_edges(extend_arrows ? B - global_scale*BIS2_combed : B, B + global_scale*BIS2_combed ,Eigen::RowVector3d(0,0,1));
  }

  if (key == '4')
  {
    // Singularities and cuts
    viewer.data().set_mesh(V, F);

    // Plot cuts
    int l_count = Seams.sum();
    Eigen::MatrixXd P1(l_count,3);
    Eigen::MatrixXd P2(l_count,3);

    for (unsigned i=0; i<Seams.rows(); ++i)
    {
      for (unsigned j=0; j<Seams.cols(); ++j)
      {
        if (Seams(i,j) != 0)
        {
          P1.row(l_count-1) = V.row(F(i,j));
          P2.row(l_count-1) = V.row(F(i,(j+1)%3));
          l_count--;
        }
      }
    }

    viewer.data().add_edges(P1, P2, Eigen::RowVector3d(1, 0, 0));

    // Plot the singularities as colored dots (red for negative, blue for positive)
    for (unsigned i=0; i<singularityIndex.size();++i)
    {
      if (singularityIndex(i) < 2 && singularityIndex(i) > 0)
        viewer.data().add_points(V.row(i),Eigen::RowVector3d(1,0,0));
      else if (singularityIndex(i) > 2)
        viewer.data().add_points(V.row(i),Eigen::RowVector3d(0,1,0));
    }

  }

  if (key == '5')
  {
    // Singularities and cuts, original field
    // Singularities and cuts
    viewer.data().set_mesh(V, F);
    viewer.data().add_edges(extend_arrows ? B - global_scale*X1_combed : B, B + global_scale*X1_combed ,Eigen::RowVector3d(1,0,0));
    viewer.data().add_edges(extend_arrows ? B - global_scale*X2_combed : B, B + global_scale*X2_combed ,Eigen::RowVector3d(0,0,1));

    // Plot cuts
    int l_count = Seams.sum();
    Eigen::MatrixXd P1(l_count,3);
    Eigen::MatrixXd P2(l_count,3);

    for (unsigned i=0; i<Seams.rows(); ++i)
    {
      for (unsigned j=0; j<Seams.cols(); ++j)
      {
        if (Seams(i,j) != 0)
        {
          P1.row(l_count-1) = V.row(F(i,j));
          P2.row(l_count-1) = V.row(F(i,(j+1)%3));
          l_count--;
        }
      }
    }

    viewer.data().add_edges(P1, P2, Eigen::RowVector3d(1, 0, 0));

    // Plot the singularities as colored dots (red for negative, blue for positive)
    for (unsigned i=0; i<singularityIndex.size();++i)
    {
      if (singularityIndex(i) < 2 && singularityIndex(i) > 0)
        viewer.data().add_points(V.row(i),Eigen::RowVector3d(1,0,0));
      else if (singularityIndex(i) > 2)
        viewer.data().add_points(V.row(i),Eigen::RowVector3d(0,1,0));
    }
  }

  if (key == '6')
  {
    // Global parametrization UV
    viewer.data().set_mesh(UV, FUV);
    viewer.data().set_uv(UV);
    viewer.data().show_lines = true;
  }

  if (key == '7')
  {
    // Global parametrization in 3D
    viewer.data().set_mesh(V, F);
    viewer.data().set_uv(UV,FUV);
    viewer.data().show_texture = true;
  }

  if (key == '8')
  {
    // Global parametrization in 3D with seams
    viewer.data().set_mesh(V, F);
    viewer.data().set_uv(UV_seams,FUV_seams);
    viewer.data().show_texture = true;
  }

  viewer.data().set_colors(Eigen::RowVector3d(1,1,1));

  // Replace the standard texture with an integer shift invariant texture
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> texture_R, texture_G, texture_B;
  line_texture(texture_R, texture_G, texture_B);
  viewer.data().set_texture(texture_R, texture_B, texture_G);

  viewer.core.align_camera_center(viewer.data().V,viewer.data().F);

  return false;
}

int main(int argc, char *argv[])
{
  using namespace Eigen;

  // Load a mesh in OFF format
  igl::readOFF(TUTORIAL_SHARED_PATH "/3holes.off", V, F);

  double gradient_size = 50;
  double iter = 0;
  double stiffness = 5.0;
  bool direct_round = 0;

  // Compute face barycenters
  igl::barycenter(V, F, B);

  // Compute scale for visualizing fields
  global_scale =  .5*igl::avg_edge_length(V, F);

    // Contrain one face
    VectorXi b(1);
    b << 0;
    MatrixXd bc(1, 3);
    bc << 1, 0, 0;

    // Create a smooth 4-RoSy field
    VectorXd S;
    igl::copyleft::comiso::nrosy(V, F, b, bc, VectorXi(), VectorXd(), MatrixXd(), 4, 0.5, X1, S);

    // Find the orthogonal vector
    MatrixXd B1, B2, B3;
    igl::local_basis(V, F, B1, B2, B3);
    X2 = igl::rotate_vectors(X1, VectorXd::Constant(1, igl::PI / 2), B1, B2);

    // Always work on the bisectors, it is more general
    igl::compute_frame_field_bisectors(V, F, X1, X2, BIS1, BIS2);

    // Comb the field, implicitly defining the seams
    igl::comb_cross_field(V, F, BIS1, BIS2, BIS1_combed, BIS2_combed);

    // Find the integer mismatches
    igl::cross_field_missmatch(V, F, BIS1_combed, BIS2_combed, true, MMatch);

    // Find the singularities
    igl::find_cross_field_singularities(V, F, MMatch, isSingularity, singularityIndex);

    // Cut the mesh, duplicating all vertices on the seams
    igl::cut_mesh_from_singularities(V, F, MMatch, Seams);

    // Comb the frame-field accordingly
    igl::comb_frame_field(V, F, X1, X2, BIS1_combed, BIS2_combed, X1_combed, X2_combed);

  // Global parametrization
  igl::copyleft::comiso::miq(V,
           F,
           X1_combed,
           X2_combed,
           MMatch,
           isSingularity,
           Seams,
           UV,
           FUV,
           gradient_size,
           stiffness,
           direct_round,
           iter,
           5,
           true);

// Global parametrization (with seams, only for demonstration)
igl::copyleft::comiso::miq(V,
         F,
         X1_combed,
         X2_combed,
         MMatch,
         isSingularity,
         Seams,
         UV_seams,
         FUV_seams,
         gradient_size,
         stiffness,
         direct_round,
         iter,
         5,
         false);

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;

  // Plot the original mesh with a texture parametrization
  key_down(viewer,'7',0);

  // Launch the viewer
  viewer.callback_key_down = &key_down;
  viewer.launch();
}

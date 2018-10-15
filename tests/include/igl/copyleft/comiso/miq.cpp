#include <gtest/gtest.h>
#include <test_common.h>
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
#include <igl/PI.h>
#include <igl/serialize.h>
#include <sstream>
#include <igl/writeDMAT.h>

TEST(miq, 3_holes)
{
using namespace Eigen;

// Input mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;

// Face barycenters
Eigen::MatrixXd B;

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

// Global parametrization
Eigen::MatrixXd UV;
Eigen::MatrixXi FUV;

// Global parametrization (reference)
Eigen::MatrixXd UV_ref;
Eigen::MatrixXi FUV_ref;

// Load a mesh in OFF format
igl::readOFF(test_common::data_path("3holes.off"), V, F);

double gradient_size = 50;
double iter = 0;
double stiffness = 5.0;
bool direct_round = 0;

// Compute face barycenters
igl::barycenter(V, F, B);

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

  // Refresh the test data
  // igl::writeDMAT(test_common::data_path("3holes-miq-UV.dmat"),UV);
  // igl::writeDMAT(test_common::data_path("3holes-miq-FUV.dmat"),FUV);

  igl::readDMAT(test_common::data_path("3holes-miq-UV.dmat"),UV_ref);
  igl::readDMAT(test_common::data_path("3holes-miq-FUV.dmat"),FUV_ref);

  ASSERT_LT((UV-UV_ref).array().abs().maxCoeff() ,1e-6);
  ASSERT_LT((FUV-FUV_ref).array().abs().maxCoeff() ,1e-6);
}

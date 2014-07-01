#undef IGL_STATIC_LIBRARY
#include <igl/readOBJ.h>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
#include <igl/viewer/Viewer.h>
#include <igl/barycenter.h>
#include <igl/avg_edge_length.h>
#include <vector>
#include <igl/n_polyvector.h>
#include <igl/conjugate_frame_fields.h>
#include <stdlib.h>
#include <igl/readOFF.h>
#include <igl/jet.h>
#include <igl/quad_planarity.h>
#include <igl/planarize_quad_mesh.h>

// Input mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;

// Face barycenters
Eigen::MatrixXd B;


// Quad mesh generated from smooth field
Eigen::MatrixXd VQS;
Eigen::MatrixXi FQS;
Eigen::MatrixXi FQStri;
Eigen::MatrixXd PQS0, PQS1, PQS2, PQS3;

// Quad mesh generated from conjugate field
Eigen::MatrixXd VQC;
Eigen::MatrixXi FQC;
Eigen::MatrixXi FQCtri;
Eigen::MatrixXd VQCplan;
Eigen::MatrixXd PQC0, PQC1, PQC2, PQC3;
Eigen::MatrixXd PQCp0, PQCp1, PQCp2, PQCp3;

// Scale for visualizing the fields
double global_scale;

// Input constraints
Eigen::VectorXi b;
Eigen::MatrixXd bc;

Eigen::MatrixXd smooth_pvf;
Eigen::MatrixXd conjugate_pvf;

igl::ConjugateFFSolverData<Eigen::MatrixXd, Eigen::MatrixXi> *csdata;

bool key_down(igl::Viewer& viewer, unsigned char key, int modifier)
{
  using namespace std;
  using namespace Eigen;

  if (key <'1' || key >'3')
    return false;

  viewer.data.lines.resize(0,9);
  // Highlight in red the constrained faces
  MatrixXd C = MatrixXd::Constant(F.rows(),3,1);
  for (unsigned i=0; i<b.size();++i)
      C.row(b(i)) << 1, 0, 0;
  viewer.set_colors(C);

  if (key == '1')
  {
    // Frame field constraints
    MatrixXd F1_t = MatrixXd::Zero(F.rows(),3);
    MatrixXd F2_t = MatrixXd::Zero(F.rows(),3);

    for (unsigned i=0; i<b.size();++i)
    {
      F1_t.row(b(i)) = bc.block(i,0,1,3);
      F2_t.row(b(i)) = bc.block(i,3,1,3);
    }

    viewer.add_edges (B - global_scale*F1_t, B + global_scale*F1_t , Eigen::RowVector3d(0,0,1));
    viewer.add_edges (B - global_scale*F2_t, B + global_scale*F2_t , Eigen::RowVector3d(0,0,1));
  }

  if (key == '2')
  {
    // Interpolated result
    viewer.add_edges (B - global_scale*smooth_pvf.block(0,0,F.rows(),3),
                      B + global_scale*smooth_pvf.block(0,0,F.rows(),3),
                      Eigen::RowVector3d(0,0,1));
    viewer.add_edges (B - global_scale*smooth_pvf.block(0,3,F.rows(),3),
                      B + global_scale*smooth_pvf.block(0,3,F.rows(),3),
                      Eigen::RowVector3d(0,0,1));
  }

  if (key == '3')
  {
    // Conjugate field
    viewer.add_edges (B - global_scale*conjugate_pvf.block(0,0,F.rows(),3),
                      B + global_scale*conjugate_pvf.block(0,0,F.rows(),3),
                      Eigen::RowVector3d(0,0,1));
    viewer.add_edges (B - global_scale*conjugate_pvf.block(0,3,F.rows(),3),
                      B + global_scale*conjugate_pvf.block(0,3,F.rows(),3),
                      Eigen::RowVector3d(0,0,1));
  }

  return false;
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;

  // Load a mesh in OBJ format
  igl::readOBJ("../shared/inspired_mesh.obj", V, F);

  // Compute face barycenters
  igl::barycenter(V, F, B);

  // Compute scale for visualizing fields
  global_scale =  .4*igl::avg_edge_length(V, F);

  // Load constraints
  igl::readDMAT("../shared/inspired_mesh_b.dmat",b);
  igl::readDMAT("../shared/inspired_mesh_bc.dmat",bc);

  // Interpolate to get a smooth field
  igl::n_polyvector(V, F, b, bc, smooth_pvf);

  // Initialize conjugate field with smooth field
  csdata = new igl::ConjugateFFSolverData<Eigen::MatrixXd,Eigen::MatrixXi>(V,F);
  conjugate_pvf = smooth_pvf;

  // Optimize the field
  int conjIter = 20;
  double lambdaOrtho = .1;
  double lambdaInit = 100;
  double lambdaMultFactor = 1.01;
  bool doHardConstraints = true;
  double lambdaOut;
  VectorXi isConstrained = VectorXi::Constant(F.rows(),0);
  for (unsigned i=0; i<b.size(); ++i)
    isConstrained(b(i)) = 1;
  
  igl::conjugate_frame_fields(*csdata, isConstrained, conjugate_pvf, conjugate_pvf, conjIter, lambdaOrtho, lambdaInit, lambdaMultFactor, doHardConstraints,
                              &lambdaOut);

  // Launch the viewer
  igl::Viewer viewer;
  viewer.core.invert_normals = true;
  viewer.core.show_lines = false;
  viewer.core.show_texture = false;
  viewer.set_mesh(V, F);
  viewer.callback_key_down = &key_down;
  key_down(viewer,'3',0);
  viewer.launch();
}

#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/conjugate_frame_fields.h>
#include <igl/ConjugateFFSolverData.h>
#include <igl/dot_row.h>
#include <igl/jet.h>
#include <igl/local_basis.h>
#include <igl/n_polyvector.h>
#include <igl/readDMAT.h>
#include <igl/readOBJ.h>
#include <igl/viewer/Viewer.h>
#include <vector>
#include <cstdlib>

#include "tutorial_shared_path.h"

// Input mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;

// Face barycenters
Eigen::MatrixXd B;

// Scale for visualizing the fields
double global_scale;

// Input constraints
Eigen::VectorXi b;
Eigen::MatrixXd bc;

Eigen::MatrixXd smooth_pvf;
Eigen::MatrixXd conjugate_pvf;
Eigen::VectorXd conjugacy_s;
Eigen::VectorXd conjugacy_c;


bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
  using namespace std;
  using namespace Eigen;

  if (key <'1' || key >'5')
    return false;

  viewer.data.lines.resize(0,9);
  // Highlight in red the constrained faces
  MatrixXd C = MatrixXd::Constant(F.rows(),3,1);
  for (unsigned i=0; i<b.size();++i)
      C.row(b(i)) << 1, 0, 0;

  double maxC = std::max(conjugacy_c.maxCoeff(), conjugacy_s.maxCoeff());
  double minC = std::min(conjugacy_c.minCoeff(), conjugacy_s.minCoeff());

  Eigen::VectorXd valS = conjugacy_s;
  // Eigen::VectorXd valS = (valS.array() - minC)/(maxC-minC);
  // valS = 1 - valS.array();
  Eigen::VectorXd valC = conjugacy_c;
  // Eigen::VectorXd valC = (valC.array() - minC)/(maxC-minC);
  // valC = 1 - valC.array();
  MatrixXd CS, CC;
  igl::jet(valS, 0, 0.004, CS);
  igl::jet(valC, 0, 0.004, CC);

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

    viewer.data.add_edges(B - global_scale*F1_t, B + global_scale*F1_t , Eigen::RowVector3d(0,0,1));
    viewer.data.add_edges(B - global_scale*F2_t, B + global_scale*F2_t , Eigen::RowVector3d(0,0,1));
    viewer.data.set_colors(C);
  }

  if (key == '2')
  {
    // Interpolated result
    viewer.data.add_edges(B - global_scale*smooth_pvf.block(0,0,F.rows(),3),
                      B + global_scale*smooth_pvf.block(0,0,F.rows(),3),
                      Eigen::RowVector3d(0,0,1));
    viewer.data.add_edges(B - global_scale*smooth_pvf.block(0,3,F.rows(),3),
                      B + global_scale*smooth_pvf.block(0,3,F.rows(),3),
                      Eigen::RowVector3d(0,0,1));
    viewer.data.set_colors(C);
  }

  if (key == '3')
  {
    // Interpolated result
    viewer.data.set_colors(CS);
  }

  if (key == '4')
  {
    // Conjugate field
    viewer.data.add_edges(B - global_scale*conjugate_pvf.block(0,0,F.rows(),3),
                      B + global_scale*conjugate_pvf.block(0,0,F.rows(),3),
                      Eigen::RowVector3d(0,0,1));
    viewer.data.add_edges(B - global_scale*conjugate_pvf.block(0,3,F.rows(),3),
                      B + global_scale*conjugate_pvf.block(0,3,F.rows(),3),
                      Eigen::RowVector3d(0,0,1));
    viewer.data.set_colors(C);
  }
  if (key == '5')
  {
    // Conjugate field
    viewer.data.set_colors(CC);
  }

  return false;
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;

  // Load a mesh in OBJ format
  igl::readOBJ(TUTORIAL_SHARED_PATH "/inspired_mesh.obj", V, F);
  // Compute face barycenters
  igl::barycenter(V, F, B);

  // Local bases (needed for conjugacy)
  Eigen::MatrixXd B1, B2, B3;
  igl::local_basis(V, F, B1, B2, B3);

  // Compute scale for visualizing fields
  global_scale =  .4*igl::avg_edge_length(V, F);

  // Load constraints
  igl::readDMAT(TUTORIAL_SHARED_PATH "/inspired_mesh_b.dmat",b);
  igl::readDMAT(TUTORIAL_SHARED_PATH "/inspired_mesh_bc.dmat",bc);

  // Interpolate to get a smooth field
  igl::n_polyvector(V, F, b, bc, smooth_pvf);

  // Initialize conjugate field with smooth field
  igl::ConjugateFFSolverData<Eigen::MatrixXd, Eigen::MatrixXi> csdata(V,F);
  conjugate_pvf = smooth_pvf;


  // Optimize the field
  int conjIter = 20;
  double lambdaOrtho = .1;
  double lambdaInit = 100;
  double lambdaMultFactor = 1.01;
  bool doHardConstraints = true;
  VectorXi isConstrained = VectorXi::Constant(F.rows(),0);
  for (unsigned i=0; i<b.size(); ++i)
    isConstrained(b(i)) = 1;

  double lambdaOut = igl::conjugate_frame_fields(csdata, isConstrained, conjugate_pvf, conjugate_pvf, conjIter, lambdaOrtho, lambdaInit, lambdaMultFactor, doHardConstraints);

  // local representations of field vectors
  Eigen::Matrix<double, Eigen::Dynamic, 2> pvU, pvV;
  pvU.resize(F.rows(),2); pvV.resize(F.rows(),2);
  //smooth
  const Eigen::MatrixXd &Us = smooth_pvf.leftCols(3);
  const Eigen::MatrixXd &Vs = smooth_pvf.rightCols(3);
  pvU << igl::dot_row(Us,B1), igl::dot_row(Us,B2);
  pvV << igl::dot_row(Vs,B1), igl::dot_row(Vs,B2);
  csdata.evaluateConjugacy(pvU, pvV, conjugacy_s);
  //conjugate
  const Eigen::MatrixXd &Uc = conjugate_pvf.leftCols(3);
  const Eigen::MatrixXd &Vc = conjugate_pvf.rightCols(3);
  pvU << igl::dot_row(Uc,B1), igl::dot_row(Uc,B2);
  pvV << igl::dot_row(Vc,B1), igl::dot_row(Vc,B2);
  csdata.evaluateConjugacy(pvU, pvV, conjugacy_c);
  // Launch the viewer
  igl::viewer::Viewer viewer;
  viewer.core.invert_normals = true;
  viewer.core.show_lines = false;
  viewer.core.show_texture = false;
  viewer.data.set_mesh(V, F);
  viewer.callback_key_down = &key_down;
  key_down(viewer,'1',0);
  viewer.launch();
}

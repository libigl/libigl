#undef IGL_STATIC_LIBRARY
#include <igl/readOBJ.h>
#include <igl/readDMAT.h>
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

// Quad mesh generated from conjugate field
Eigen::MatrixXd VQC;
Eigen::MatrixXi FQC;
Eigen::MatrixXi FQCtri;
Eigen::MatrixXd VQCplan;

// Scale for visualizing the fields
double global_scale;

// Input constraints
Eigen::VectorXi isConstrained;
Eigen::MatrixXd constraints;

Eigen::MatrixXd smooth_pvf;
Eigen::MatrixXd conjugate_pvf;

igl::ConjugateFFSolverData<Eigen::MatrixXd, Eigen::MatrixXi> *csdata;

int conjIter = 2;
int totalConjIter = 0;
double lambdaOrtho = .1;
double lambdaInit = 100;
double lambdaMultFactor = 1.01;
bool doHardConstraints = true;

bool key_down(igl::Viewer& viewer, unsigned char key, int modifier)
{
  using namespace std;
  using namespace Eigen;

  if (key <'1' || key >'6')
    return false;

  viewer.clear_mesh();
  viewer.options.show_lines = false;
  viewer.options.show_texture = false;
  if (key <= '3')
  {
    viewer.set_mesh(V, F);
    // Highlight in red the constrained faces
    MatrixXd C = MatrixXd::Constant(F.rows(),3,1);
    for (unsigned i=0; i<F.rows();++i)
      if (isConstrained[i])
        C.row(i) << 1, 0, 0;
    viewer.set_colors(C);
  }

  if (key == '1')
  {
    // Frame field constraints
    MatrixXd F1_t = MatrixXd::Zero(F.rows(),3);
    MatrixXd F2_t = MatrixXd::Zero(F.rows(),3);
    for (unsigned i=0; i<F.rows();++i)
      if (isConstrained[i])
      {
        F1_t.row(i) = constraints.block(i,0,1,3);
        F2_t.row(i) = constraints.block(i,3,1,3);
      }
    viewer.add_edges (B - global_scale*F1_t, B + global_scale*F1_t , Eigen::RowVector3d(0,0,1));
    viewer.add_edges (B - global_scale*F2_t, B + global_scale*F2_t , Eigen::RowVector3d(0,0,1));

  }
  if (key == '2')
  {
    viewer.add_edges (B - global_scale*smooth_pvf.block(0,0,F.rows(),3),
                      B + global_scale*smooth_pvf.block(0,0,F.rows(),3),
                      Eigen::RowVector3d(0,1,0));
    viewer.add_edges (B - global_scale*smooth_pvf.block(0,3,F.rows(),3),
                      B + global_scale*smooth_pvf.block(0,3,F.rows(),3),
                      Eigen::RowVector3d(0,1,0));
  }

  if (key == '3')
  {
    if (totalConjIter <50)
    {
      double lambdaOut;
      igl::conjugate_frame_fields(*csdata, isConstrained, conjugate_pvf, conjugate_pvf, conjIter, lambdaOrtho, lambdaInit, lambdaMultFactor, doHardConstraints,
                                  &lambdaOut);
      totalConjIter += 2;
      lambdaInit = lambdaOut;
    }
    viewer.add_edges (B - global_scale*conjugate_pvf.block(0,0,F.rows(),3),
                      B + global_scale*conjugate_pvf.block(0,0,F.rows(),3),
                      Eigen::RowVector3d(0,1,0));
    viewer.add_edges (B - global_scale*conjugate_pvf.block(0,3,F.rows(),3),
                      B + global_scale*conjugate_pvf.block(0,3,F.rows(),3),
                      Eigen::RowVector3d(0,1,0));
  }

  if (key == '4')
  {
    viewer.set_mesh(VQS, FQStri);
    
    // show planarity
    VectorXd planarity;
    igl::quad_planarity( VQS, FQS, planarity);
    MatrixXd Ct;
    igl::jet(planarity, 0, 0.02, Ct);
    MatrixXd C(FQStri.rows(),3);
    C << Ct, Ct;
    viewer.set_colors(C);
  }
  
  if (key == '5')
  {
    viewer.set_mesh(VQC, FQCtri);
    
    // show planarity
    VectorXd planarity;
    igl::quad_planarity( VQC, FQC, planarity);
    MatrixXd Ct;
    igl::jet(planarity, 0, 0.02, Ct);
    MatrixXd C(FQCtri.rows(),3);
    C << Ct, Ct;
    viewer.set_colors(C);
  }
  
  if (key == '6')
  {
    igl ::planarize_quad_mesh(VQC, FQC, 50, 0.01, VQCplan);
    viewer.set_mesh(VQCplan, FQCtri);

    // show planarity
    VectorXd planarity;
    igl::quad_planarity( VQCplan, FQC, planarity);
    MatrixXd Ct;
    igl::jet(planarity, 0, 0.02, Ct);
    MatrixXd C(FQCtri.rows(),3);
    C << Ct, Ct;
    viewer.set_colors(C);
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
  global_scale =  .2*igl::avg_edge_length(V, F);

  // Load constraints
  MatrixXd temp;
  igl::readDMAT("../shared/inspired_mesh.dmat",temp);
  isConstrained = temp.block(0,0,temp.rows(),1).cast<int>();
  constraints = temp.block(0,1,temp.rows(),temp.cols()-1);

  // Interpolate to get a smooth field
  igl::n_polyvector(V, F, isConstrained, constraints, smooth_pvf);

  // Initialize conjugate field with smooth field
  csdata = new igl::ConjugateFFSolverData<Eigen::MatrixXd,Eigen::MatrixXi>(V,F);
  conjugate_pvf = smooth_pvf;

  // Load quad mesh generated by smooth field
  igl::readOFF("../shared/inspired_mesh_quads_Smooth.off", VQS, FQS);
  FQStri.resize(2*FQS.rows(), 3);
  FQStri <<  FQS.col(0),FQS.col(1),FQS.col(2),
             FQS.col(2),FQS.col(3),FQS.col(0);

  // Load quad mesh generated by conjugate field
  igl::readOFF("../shared/inspired_mesh_quads_Conjugate.off", VQC, FQC);
  FQCtri.resize(2*FQC.rows(), 3);
  FQCtri <<  FQC.col(0),FQC.col(1),FQC.col(2),
             FQC.col(2),FQC.col(3),FQC.col(0);



  igl::Viewer viewer;

  // Plot the original mesh with a texture parametrization
  key_down(viewer,'1',0);

  // Launch the viewer
  viewer.callback_key_down = &key_down;
  viewer.launch();
}

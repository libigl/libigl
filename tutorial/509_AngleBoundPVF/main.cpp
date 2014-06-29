#undef IGL_STATIC_LIBRARY
#include <igl/readOBJ.h>
#include <igl/readDMAT.h>
#include <igl/viewer/Viewer.h>
#include <igl/barycenter.h>
#include <igl/avg_edge_length.h>
#include <vector>
#include <igl/n_polyvector.h>
#include <igl/angle_bound_frame_fields.h>
#include <stdlib.h>
#include <igl/jet.h>

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
Eigen::MatrixXd PQC0, PQC1, PQC2, PQC3;

// Scale for visualizing the fields
double global_scale;

// Input constraints
Eigen::VectorXi isConstrained;
Eigen::MatrixXd constraints;

Eigen::MatrixXd smooth_pvf;
Eigen::MatrixXd angle_bound_pvf;

igl::AngleBoundFFSolverData<Eigen::MatrixXd, Eigen::MatrixXi> *csdata;

int conjIter = 2;
int totalConjIter = 0;
double lambdaOrtho = .1;
double lambdaInit = 100;
double lambdaMultFactor = 1.5;
bool doHardConstraints = false;

bool showAngles = true;
int curr_key = 0;

void computeAngles(const Eigen::MatrixXd &ff, Eigen::VectorXd &angles)
{
 angles.resize(ff.rows(),1);
  int num =0;
 for (int i =0; i<ff.rows(); ++i)
 {
  Eigen::RowVector3d u = (ff.block(i,0,1,3)); u.normalize();
  Eigen::RowVector3d v = (ff.block(i,3,1,3)); v.normalize();
  double s = (u.cross(v)).norm();
  double c = fabs(u.dot(v));
  angles[i] = atan2(s,c);
   num += (angles[i]<70*M_PI/180);
  }
  std::cerr<<"out of bound:"<<num<<std::endl;
}

void getAngleColor(const Eigen::MatrixXd &ff, Eigen::MatrixXd &C)
{
  Eigen::VectorXd angles;
  computeAngles(ff, angles);
  Eigen::VectorXd val = 0.5*M_PI*Eigen::VectorXd::Ones(angles.rows(),1)-angles;
  igl::jet(val, 0, 20*M_PI/180., C);
}

bool key_down(igl::Viewer& viewer, unsigned char key, int modifier)
{
  using namespace std;
  using namespace Eigen;

  // Highlight in red the constrained faces
  MatrixXd CC = MatrixXd::Constant(F.rows(),3,1);
  for (unsigned i=0; i<F.rows();++i)
    if (isConstrained[i])
      CC.row(i) << 1, 0, 0;

  if (key == 'c' || key == 'C')
  {
    showAngles = !showAngles;
    if (curr_key == 2)
    {
      MatrixXd C = CC;
      if (showAngles)
        getAngleColor(smooth_pvf, C);
      viewer.set_colors(C);
    }
    else if (curr_key == 3)
    {
      MatrixXd C = CC;
      if (showAngles)
        getAngleColor(angle_bound_pvf, C);
      viewer.set_colors(C);
    }
    return false;
  }
  
  if (key <'1' || key >'5')
  {
    return false;
  }

  viewer.clear_mesh();
  viewer.options.show_lines = false;
  viewer.options.show_texture = false;

  if (key == '1')
  {
    viewer.set_mesh(V, F);
    viewer.set_colors(CC);
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
    curr_key = 1;
  }
  if (key == '2')
  {
    viewer.set_mesh(V, F);
    viewer.add_edges (B - global_scale*smooth_pvf.block(0,0,F.rows(),3),
                      B + global_scale*smooth_pvf.block(0,0,F.rows(),3),
                      Eigen::RowVector3d(0,1,0));
    viewer.add_edges (B - global_scale*smooth_pvf.block(0,3,F.rows(),3),
                      B + global_scale*smooth_pvf.block(0,3,F.rows(),3),
                      Eigen::RowVector3d(0,1,0));
    MatrixXd C = CC;
    if (showAngles)
      getAngleColor(smooth_pvf, C);
    viewer.set_colors(C);
    curr_key = 2;
  }

  if (key == '3')
  {
    viewer.set_mesh(V, F);
    if (totalConjIter <50)
    {
      double lambdaOut;
      igl::angle_bound_frame_fields(*csdata,
                                     70,
                                     isConstrained,
                                     angle_bound_pvf,
                                     angle_bound_pvf,
                                     conjIter,
                                     lambdaInit,
                                     lambdaMultFactor,
                                     doHardConstraints,
                                     &lambdaOut);
      totalConjIter += conjIter;
      lambdaInit = lambdaOut;
    }
    viewer.add_edges (B - global_scale*angle_bound_pvf.block(0,0,F.rows(),3),
                      B + global_scale*angle_bound_pvf.block(0,0,F.rows(),3),
                      Eigen::RowVector3d(0,1,0));
    viewer.add_edges (B - global_scale*angle_bound_pvf.block(0,3,F.rows(),3),
                      B + global_scale*angle_bound_pvf.block(0,3,F.rows(),3),
                      Eigen::RowVector3d(0,1,0));
    MatrixXd C = CC;
    if (showAngles)
      getAngleColor(angle_bound_pvf, C);
    viewer.set_colors(C);
    curr_key = 3;
  }

  if (key == '4')
  {
    viewer.set_mesh(VQS, FQStri);
    viewer.add_edges (PQS0, PQS1, Eigen::RowVector3d(0,0,0));
    viewer.add_edges (PQS1, PQS2, Eigen::RowVector3d(0,0,0));
    viewer.add_edges (PQS2, PQS3, Eigen::RowVector3d(0,0,0));
    viewer.add_edges (PQS3, PQS0, Eigen::RowVector3d(0,0,0));
    curr_key = 4;
  }

  if (key == '5')
  {
    viewer.set_mesh(VQC, FQCtri);
    viewer.add_edges (PQC0, PQC1, Eigen::RowVector3d(0,0,0));
    viewer.add_edges (PQC1, PQC2, Eigen::RowVector3d(0,0,0));
    viewer.add_edges (PQC2, PQC3, Eigen::RowVector3d(0,0,0));
    viewer.add_edges (PQC3, PQC0, Eigen::RowVector3d(0,0,0));
    curr_key = 5;
  }


  return false;
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  // Load a mesh in OBJ format
  igl::readOBJ("../shared/teddy.obj", V, F);

  // Compute face barycenters
  igl::barycenter(V, F, B);

  // Compute scale for visualizing fields
  global_scale =  .2*igl::avg_edge_length(V, F);

  // Load constraints
  MatrixXd temp;
  igl::readDMAT("../shared/teddy.dmat",temp);
  isConstrained = temp.block(0,0,temp.rows(),1).cast<int>();
  constraints = temp.block(0,1,temp.rows(),temp.cols()-1);

  // Interpolate to get a smooth field
  igl::n_polyvector(V, F, isConstrained, constraints, smooth_pvf);

  // Initialize conjugate field with smooth field
  csdata = new igl::AngleBoundFFSolverData<Eigen::MatrixXd,Eigen::MatrixXi>(V,F);
  angle_bound_pvf = smooth_pvf;

  // Load quad mesh generated by smooth field
  igl::readOBJ("../shared/teddy_smooth_remeshed.obj", VQS, FQS);
  FQStri.resize(2*FQS.rows(), 3);
  FQStri <<  FQS.col(0),FQS.col(1),FQS.col(2),
             FQS.col(2),FQS.col(3),FQS.col(0);

  // Load quad mesh generated by conjugate field
  igl::readOBJ("../shared/teddy_angle_bound_remeshed.obj", VQC, FQC);
  FQCtri.resize(2*FQC.rows(), 3);
  FQCtri <<  FQC.col(0),FQC.col(1),FQC.col(2),
             FQC.col(2),FQC.col(3),FQC.col(0);

  igl::slice( VQS, FQS.col(0), 1, PQS0);
  igl::slice( VQS, FQS.col(1), 1, PQS1);
  igl::slice( VQS, FQS.col(2), 1, PQS2);
  igl::slice( VQS, FQS.col(3), 1, PQS3);

  igl::slice( VQC, FQC.col(0), 1, PQC0);
  igl::slice( VQC, FQC.col(1), 1, PQC1);
  igl::slice( VQC, FQC.col(2), 1, PQC2);
  igl::slice( VQC, FQC.col(3), 1, PQC3);

  igl::Viewer viewer;

  // Plot the original mesh with a texture parametrization
  key_down(viewer,'1',0);

  // Launch the viewer
  viewer.callback_key_down = &key_down;
  viewer.launch();
}

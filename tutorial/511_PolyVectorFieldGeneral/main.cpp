#include <igl/readOBJ.h>
#include <igl/readDMAT.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/barycenter.h>
#include <igl/avg_edge_length.h>
#include <vector>
#include <igl/n_polyvector_general.h>
#include <igl/n_polyvector.h>
#include <igl/local_basis.h>
#include <igl/writeOFF.h>
#include <stdlib.h>
#include <igl/jet.h>
#include <fstream>
#include <iostream>
// Input mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;

// Per face bases
Eigen::MatrixXd B1,B2,B3;

// Face barycenters
Eigen::MatrixXd B;

// Scale for visualizing the fields
double global_scale;

// Random length factor
double rand_factor = 5;

Eigen::VectorXi samples;

void readSamples(const std::string &fname, Eigen::VectorXi &samples)
{
  int numSamples;
  FILE *fp = fopen(fname.c_str(),"r");
  if (fscanf(fp, "%d", &numSamples)!=1)
  {
    fclose(fp);
    return;
  }
  samples.resize(numSamples,1);
  int vali;
  for (int i =0; i<numSamples; ++i)
  {
    if (fscanf(fp, "%d", &vali)!=1 || vali<0)
    {
      fclose(fp);
      samples.resize(0,1);
      return;
    }
    samples[i]=vali;
  }
  fclose(fp);
  
}

// Create a random set of tangent vectors
Eigen::VectorXd random_constraints(const
                                   Eigen::VectorXd& b1, const
                                   Eigen::VectorXd& b2, int n)
{
  Eigen::VectorXd r(n*3);
  for (unsigned i=0; i<n;++i)
  {
    double a = (double(rand())/RAND_MAX)*2*M_PI;
    double s = 1 + ((double(rand())/RAND_MAX)) * rand_factor;
    Eigen::Vector3d t = s * (cos(a) * b1 + sin(a) * b2);
    r.block(i*3,0,3,1) = t;
  }
  return r;
}

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  using namespace std;
  using namespace Eigen;

  if (key <'1' || key >'9')
    return false;

  viewer.data().lines.resize(0,9);

  int num = key  - '0';

  // Interpolate
  std::cerr << "Interpolating " << num  << "-PolyVector field" << std::endl;

  VectorXi b(3);
  b << 1511, 603, 506;

  int numConstraintsToGenerate;
  // if it's not a 2-PV or a 1-PV, include a line direction (2 opposite vectors)
  // in the field
  if (num>=5)
    numConstraintsToGenerate  = num-2;
  else
    if (num>=3)
      numConstraintsToGenerate  = num-1;
    else
      numConstraintsToGenerate  = num;


  MatrixXd bc(b.size(),numConstraintsToGenerate*3);
  for (unsigned i=0; i<b.size(); ++i)
  {
    VectorXd t = random_constraints(B1.row(b(i)),B2.row(b(i)),numConstraintsToGenerate);
    bc.row(i) = t;
  }
  VectorXi rootsIndex(num);
  for (int i =0; i<numConstraintsToGenerate; ++i)
    rootsIndex[i] = i+1;
  if (num>=5)
    rootsIndex[num-2] = -2;
    if (num>=3)
      rootsIndex[num-1] = -1;

  // Interpolated PolyVector field
  Eigen::MatrixXd pvf;
  igl::n_polyvector_general(V, F, b, bc, rootsIndex, pvf);

  ofstream ofs;
  ofs.open("pvf.txt", ofstream::out);
  ofs<<pvf;
  ofs.close();
  igl::writeOFF("pvf.off",V,F);
  
  // Highlight in red the constrained faces
  MatrixXd C = MatrixXd::Constant(F.rows(),3,1);
  for (unsigned i=0; i<b.size();++i)
    C.row(b(i)) << 1, 0, 0;
  viewer.data().set_colors(C);

  for (int n=0; n<num; ++n)
  {
//    const MatrixXd &VF = pvf.block(0,n*3,F.rows(),3);
    MatrixXd VF = MatrixXd::Zero(F.rows(),3);
    for (unsigned i=0; i<b.size(); ++i)
      VF.row(b[i]) = pvf.block(b[i],n*3,1,3);
    
    for (int i=0; i<samples.rows(); ++i)
      VF.row(samples[i]) = pvf.block(samples[i],n*3,1,3);
    
    VectorXd c = VF.rowwise().norm();
    MatrixXd C2;
    igl::jet(c,1,1+rand_factor,C2);
    // viewer.data().add_edges(B - global_scale*VF, B + global_scale*VF , C2);
    viewer.data().add_edges(B, B + global_scale*VF , C2);
  }


  return false;
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  // Load a mesh in OBJ format
  igl::readOBJ(TUTORIAL_SHARED_PATH "/snail.obj", V, F);
  readSamples(TUTORIAL_SHARED_PATH "/snail.samples.0.2", samples);

  // Compute local basis for faces
  igl::local_basis(V,F,B1,B2,B3);

  // Compute face barycenters
  igl::barycenter(V, F, B);

  // Compute scale for visualizing fields
  global_scale =  .1*igl::avg_edge_length(V, F);

  // Make the example deterministic
  srand(0);

  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.callback_key_down = &key_down;
  viewer.data().show_lines = false;
  key_down(viewer,'3',0);
  std::cerr << " *** Press keys 1-9 to select number of vectors per point. ***" << std::endl;
  
  viewer.launch();
}

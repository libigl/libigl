#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/jet.h>
#include <igl/local_basis.h>
#include <igl/n_polyvector.h>
#include <igl/readDMAT.h>
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "tutorial_shared_path.h"

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

  if (key <'1' || key >'8')
    return false;

  viewer.data().lines.resize(0,9);

  int num = key  - '0';

  // Interpolate
  std::cerr << "Interpolating " << num * 2 << "-PolyVector field" << std::endl;

  VectorXi b(4);
  b << 4550, 2321, 5413, 5350;

  MatrixXd bc(b.size(),num*3);
  for (unsigned i=0; i<b.size(); ++i)
  {
    VectorXd t = random_constraints(B1.row(b(i)),B2.row(b(i)),num);
    bc.row(i) = t;
  }

  // Interpolated PolyVector field
  Eigen::MatrixXd pvf;
  igl::n_polyvector(V, F, b, bc, pvf);

  // Highlight in red the constrained faces
  MatrixXd C = MatrixXd::Constant(F.rows(),3,1);
  for (unsigned i=0; i<b.size();++i)
    C.row(b(i)) << 1, 0, 0;
  viewer.data().set_colors(C);

  for (int n=0; n<num; ++n)
  {
    MatrixXd VF = MatrixXd::Zero(F.rows(),3);
    for (unsigned i=0; i<b.size(); ++i)
      VF.row(b[i]) = bc.row(i);

    for (int i=0; i<samples.rows(); ++i)
      VF.row(samples[i]) = pvf.block(samples[i],n*3,1,3);
    // MatrixXd VF = pvf.block(0,n*3,F.rows(),3);

    VectorXd c = VF.rowwise().norm();
    MatrixXd C2;
    igl::jet(c,1,1+rand_factor,C2);
    viewer.data().add_edges(B - global_scale*VF, B + global_scale*VF , C2);
  }

  return false;
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  // Load a mesh in OBJ format
  igl::readOBJ(TUTORIAL_SHARED_PATH "/lilium.obj", V, F);
  readSamples(TUTORIAL_SHARED_PATH "/lilium.samples.0.2", samples);

  // Compute local basis for faces
  igl::local_basis(V,F,B1,B2,B3);

  // Compute face barycenters
  igl::barycenter(V, F, B);

  // Compute scale for visualizing fields
  global_scale =  .2*igl::avg_edge_length(V, F);

  // Make the example deterministic
  srand(0);

  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.callback_key_down = &key_down;
  viewer.data().show_lines = false;
  viewer.data().line_width = 10000;// this does not work, why?
  key_down(viewer,'2',0);


  viewer.launch();
}

#include "parse_rhs.h"

#include <igl/matlab/mexStream.h>
#include <igl/matlab/MatlabWorkspace.h>
#include <igl/embree/EmbreeIntersector.h>
#include <igl/embree/ambient_occlusion.h>
#include <igl/matlab_format.h>

#include <igl/read.h>
#include <igl/per_vertex_normals.h>

#include <mex.h>

#include <iostream>
#include <string>

void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
  // This is useful for debugging whether Matlab is caching the mex binary
  //mexPrintf("%s %s\n",__TIME__,__DATE__);
  igl::mexStream mout;
  std::streambuf *outbuf = std::cout.rdbuf(&mout);

  using namespace std;
  using namespace Eigen;
  using namespace igl;

  MatrixXd V,P,N;
  VectorXd S;
  MatrixXi F;
  int num_samples;
  parse_rhs(nrhs,prhs,V,F,P,N,num_samples);
  // Prepare left-hand side
  nlhs = 1;

  //read("../shared/cheburashka.obj",V,F);
  //P = V;
  //per_vertex_normals(V,F,N);
  EmbreeIntersector<::MatrixXd::Scalar,::MatrixXi::Scalar> ei;
  ei = EmbreeIntersector<MatrixXd::Scalar,MatrixXi::Scalar>(V,F);
  ambient_occlusion(ei,P,N,num_samples,S);
  MatlabWorkspace mw;
  mw.save(V,"V");
  mw.save(P,"P");
  mw.save(N,"N");
  mw.save_index(F,"F");
  mw.save(S,"S");
  mw.write("out.mat");

  plhs[0] = mxCreateDoubleMatrix(S.rows(),S.cols(), mxREAL);
  copy(S.data(),S.data()+S.size(),mxGetPr(plhs[0]));

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
}

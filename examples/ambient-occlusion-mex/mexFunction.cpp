#include "parse_rhs.h"

#include <igl/matlab/MexStream.h>
#include <igl/embree/ambient_occlusion.h>

#include <igl/read_triangle_mesh.h>
#include <igl/per_vertex_normals.h>

#include <mex.h>

#include <iostream>
#include <string>

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
  // This is useful for debugging whether Matlab is caching the mex binary
  //mexPrintf("%s %s\n",__TIME__,__DATE__);
  igl::MexStream mout;
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

  //read_triangle_mesh("../shared/cheburashka.off",V,F);
  //P = V;
  //per_vertex_normals(V,F,N);
  ambient_occlusion(V,F,P,N,num_samples,S);
  //MatlabWorkspace mw;
  //mw.save(V,"V");
  //mw.save(P,"P");
  //mw.save(N,"N");
  //mw.save_index(F,"F");
  //mw.save(S,"S");
  //mw.write("out.mat");

  plhs[0] = mxCreateDoubleMatrix(S.rows(),S.cols(), mxREAL);
  copy(S.data(),S.data()+S.size(),mxGetPr(plhs[0]));

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
}

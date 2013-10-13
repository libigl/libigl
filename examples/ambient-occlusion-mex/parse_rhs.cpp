#ifdef MEX
#include "parse_rhs.h"
#include <algorithm>
#include <functional>

void parse_rhs(
  const int nrhs, 
  const mxArray *prhs[], 
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F,
  Eigen::MatrixXd & P,
  Eigen::MatrixXd & N,
  int & num_samples)
{
  using namespace std;
  if(nrhs < 5)
  {
    mexErrMsgTxt("nrhs < 5");
  }

  const int dim = mxGetN(prhs[0]);
  if(dim != 3)
  {
    mexErrMsgTxt("Mesh vertex list must be #V by 3 list of vertex positions");
  }
  if(dim != (int)mxGetN(prhs[1]))
  {
   mexErrMsgTxt("Mesh facet size must be 3");
  }
  if(mxGetN(prhs[2]) != dim)
  {
    mexErrMsgTxt("Point list must be #P by 3 list of origin locations");
  }
  if(mxGetN(prhs[3]) != dim)
  {
    mexErrMsgTxt("Normal list must be #P by 3 list of origin normals");
  }
  if(mxGetN(prhs[4]) != 1 || mxGetM(prhs[4]) != 1)
  {
    mexErrMsgTxt("Number of samples must be scalar.");
  }


  V.resize(mxGetM(prhs[0]),mxGetN(prhs[0]));
  copy(mxGetPr(prhs[0]),mxGetPr(prhs[0])+V.size(),V.data());
  F.resize(mxGetM(prhs[1]),mxGetN(prhs[1]));
  copy(mxGetPr(prhs[1]),mxGetPr(prhs[1])+F.size(),F.data());
  F.array() -= 1;
  P.resize(mxGetM(prhs[2]),mxGetN(prhs[2]));
  copy(mxGetPr(prhs[2]),mxGetPr(prhs[2])+P.size(),P.data());
  N.resize(mxGetM(prhs[3]),mxGetN(prhs[3]));
  copy(mxGetPr(prhs[3]),mxGetPr(prhs[3])+N.size(),N.data());
  if(*mxGetPr(prhs[4]) != (int)*mxGetPr(prhs[4]))
  {
    mexErrMsgTxt("Number of samples should be non negative integer.");
  }
  num_samples = (int) *mxGetPr(prhs[4]);
}
#endif

#include "mex.h"

#define IGL_HEADER_ONLY
#include <Eigen/Core>
#include <igl/readOBJ.h>


using namespace std;
using namespace Eigen;

extern void _main();

Eigen::MatrixXd readMatrix(const mxArray* mat)
{
  double* ptr = mxGetPr(mat);

  int m = mxGetM(mat);
  int n = mxGetN(mat);

  Eigen::MatrixXd	V;
  V.resize(m,n);
  for(int j=0;j<n;j++)
    for(int i=0;i<m;++i)
      V(i,j) = *(ptr++);

  return V;
}

mxArray* writeMatrix(const MatrixXd& mat)
{
  mxArray* ret = 0;

  if (mat.size() == 0)
  {
    ret =	mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
  }
  else
  {
    int M = mat.rows();
    int N = mat.cols();

    ret = mxCreateNumericMatrix(M, N, mxDOUBLE_CLASS, mxREAL);
    double* pointer = mxGetPr(ret);

    /* Copy data into the mxArray */
    int count = 0;
    for ( int j = 0; j < N; ++j )
      for (int i = 0; i < M; ++i)
        pointer[count++] = mat(i,j);
  }

  return ret;
}

void mexFunction(
     int          nlhs,
     mxArray      *plhs[],
     int          nrhs,
     const mxArray *prhs[]
     )
{
  /* Check for proper number of arguments */

  if (nrhs != 1) {
  mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin",
      "readOBJ requires 1 input arguments, the path of the file to open");
  }
  else if (nlhs != 2) {
  mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout",
      "readOBJ generates two output argument, V and F.");
  }

  // Read the file path
  char* file_path = mxArrayToString(prhs[0]);

  MatrixXd V;
  MatrixXi F;

  // Read the mesh
  igl::readOBJ(file_path,V,F);

  // Matlab is indexing matrices from 1
  F = F.array() + 1;

  // Return the matrices to matlab, after converting them to double
  plhs[0] = writeMatrix(V);
  plhs[1] = writeMatrix(F.cast<double>());

  return;
}

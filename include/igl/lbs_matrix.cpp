#include "lbs_matrix.h"

IGL_INLINE void igl::lbs_matrix(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXd & W,
  Eigen::SparseMatrix<double>& M)
{
  // number of mesh vertices
  int n = V.rows();
  assert(n == W.rows());
  // dimension of mesh
  int dim = V.cols();
  // number of handles
  int m = W.cols();

  M.resize(n*dim,m*dim*(dim+1));

  // loop over coordinates of mesh vertices
  for(int x = 0; x < dim; x++)
  {
    // loop over mesh vertices
    for(int j = 0; j < n; j++)
    {
      // loop over handles
      for(int i = 0; i < m; i++)
      {
        // loop over cols of affine transformations
        for(int c = 0; c < (dim+1); c++)
        {
          double value = W(j,i);
          if(c<dim)
          {
            value *= V(j,c);
          }
          M.insert(x*n + j,x*m + c*m*dim + i) = value;
        }
      }
    }
  }

  M.makeCompressed();
}

IGL_INLINE void igl::lbs_matrix(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXd & W,
  Eigen::MatrixXd & M)
{
  // number of mesh vertices
  int n = V.rows();
  assert(n == W.rows());
  // dimension of mesh
  int dim = V.cols();
  // number of handles
  int m = W.cols();
  M.resize(n*dim,m*dim*(dim+1));

  // loop over coordinates of mesh vertices
  for(int x = 0; x < dim; x++)
  {
    // loop over mesh vertices
    for(int j = 0; j < n; j++)
    {
      // loop over handles
      for(int i = 0; i < m; i++)
      {
        // loop over cols of affine transformations
        for(int c = 0; c < (dim+1); c++)
        {
          double value = W(j,i);
          if(c<dim)
          {
            value *= V(j,c);
          }
          M(x*n + j,x*m + c*m*dim + i) = value;
        }
      }
    }
  }
}

IGL_INLINE void igl::lbs_matrix(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXd & W,
  const Eigen::MatrixXi & WI,
  Eigen::SparseMatrix<double>& M)
{
  // number of mesh vertices
  int n = V.rows();
  assert(n == W.rows());
  assert(n == WI.rows());
  // dimension of mesh
  int dim = V.cols();
  // number of handles
  int m = WI.maxCoeff()+1;
  // max number of influencing handles
  int k = W.cols();
  assert(k == WI.cols());

  M.resize(n*dim,m*dim*(dim+1));

  // loop over coordinates of mesh vertices
  for(int x = 0; x < dim; x++)
  {
    // loop over mesh vertices
    for(int j = 0; j < n; j++)
    {
      // loop over handles
      for(int i = 0; i < k; i++)
      {
        // loop over cols of affine transformations
        for(int c = 0; c < (dim+1); c++)
        {
          double value = W(j,i);
          if(c<dim)
          {
            value *= V(j,c);
          }
          if(value != 0)
          {
            M.insert(x*n + j,x*m + c*m*dim + WI(j,i)) = value;
          }
        }
      }
    }
  }

  M.makeCompressed();
}


IGL_INLINE void igl::lbs_matrix(
  const Eigen::MatrixXd & V, 
  const Eigen::MatrixXd & W,
  const Eigen::MatrixXi & WI,
  Eigen::MatrixXd & M)
{
  // Cheapskate wrapper
  using namespace Eigen;
  SparseMatrix<double> sM;
  lbs_matrix(V,W,WI,sM);
  M = MatrixXd(sM);
}

#include "sample_edges.h"

IGL_INLINE void igl::sample_edges(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const int k,
  Eigen::MatrixXd & S)
{
  using namespace Eigen;
  // Resize output
  S.resize(V.rows() + k * E.rows(),V.cols());
  // Copy V at front of S
  S.block(0,0,V.rows(),V.cols()) = V;

  // loop over edges
  for(int i = 0;i<E.rows();i++)
  {
    VectorXd tip = V.row(E(i,0));
    VectorXd tail = V.row(E(i,1));
    for(int s=0;s<k;s++)
    {
      double f = double(s+1)/double(k+1);
      S.row(V.rows()+k*i+s) = f*tail + (1.0-f)*tip;
    }
  }
}

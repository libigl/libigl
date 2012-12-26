#include "bounding_box_diagonal.h"
#include "mat_max.h"
#include "mat_min.h"
#include <cmath>

IGL_INLINE double igl::bounding_box_diagonal(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F)
{
  using namespace igl;
  using namespace Eigen;
  VectorXd maxV,minV;
  VectorXi maxVI,minVI;
  mat_max(V,1,maxV,maxVI);
  mat_min(V,1,minV,minVI);
  return sqrt((maxV-minV).array().square().sum());
}

#ifndef IGL_HEADER_ONLY
// Explicit template instanciation
#endif

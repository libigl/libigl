#include "unique_simplices.h"
#include "sort.h"
#include "unique.h"

IGL_INLINE void igl::unique_simplices(
  const Eigen::MatrixXi & F,
  Eigen::MatrixXi & FF)
{
  using namespace Eigen;
  using namespace igl;
  // Sort each face
  MatrixXi sortF, unusedI;
  igl::sort(F,2,1,sortF,unusedI);
  // Find unique faces
  VectorXi IA,IC;
  MatrixXi C;
  igl::unique_rows(sortF,C,IA,IC);
  FF.resize(IA.size(),F.cols());
  // Copy into output
  for(int i = 0;i<IA.rows();i++)
  {
    FF.row(i) = F.row(IA(i));
  }
}

#include "ears.h"
#include "on_boundary.h"
#include "find.h"
#include "mat_min.h"
#include <cassert>

template <
  typename DerivedF,
  typename Derivedear,
  typename Derivedear_opp>
IGL_INLINE void igl::ears(
  const Eigen::MatrixBase<DerivedF> & F,
  Eigen::PlainObjectBase<Derivedear> & ear,
  Eigen::PlainObjectBase<Derivedear_opp> & ear_opp)
{
  assert(F.cols() == 3 && "F should contain triangles");
  Eigen::Array<bool,Eigen::Dynamic,3> B;
  {
    Eigen::Array<bool,Eigen::Dynamic,1> I;
    on_boundary(F,I,B);
  }
  find(B.rowwise().count() == 2,ear);
  Eigen::Array<bool,Eigen::Dynamic,3> Bear;
  slice(B,ear,1,Bear);
  Eigen::Array<bool,Eigen::Dynamic,1> M;
  mat_min(Bear,2,M,ear_opp);
}

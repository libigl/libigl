#include "directed_edge_parents.h"
#include "slice_into.h"
#include "slice.h"
#include "colon.h"
#include "setdiff.h"
#include <algorithm>

template <typename DerivedE, typename DerivedP>
IGL_INLINE void igl::directed_edge_parents(
  const Eigen::PlainObjectBase<DerivedE> & E,
  Eigen::PlainObjectBase<DerivedP> & P)
{
  using namespace Eigen;
  using namespace std;
  VectorXi I = VectorXi::Constant(E.maxCoeff()+1,1,-1);
  //I(E.col(1)) = 0:E.rows()-1
  slice_into(colon<int>(0,E.rows()-1),E.col(1).eval(),I);
  VectorXi roots,_;
  setdiff(E.col(0).eval(),E.col(1).eval(),roots,_);
  for_each(roots.data(),roots.data()+roots.size(),[&](int r){I(r)=-1;});
  slice(I,E.col(0).eval(),P);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instanciation
template void igl::directed_edge_parents<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif

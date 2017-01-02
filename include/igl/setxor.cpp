#include "setxor.h"
#include "setdiff.h"
#include "setunion.h"
#include "slice.h"

template <
  typename DerivedA,
  typename DerivedB,
  typename DerivedC,
  typename DerivedIA,
  typename DerivedIB>
IGL_INLINE void igl::setxor(
  const Eigen::DenseBase<DerivedA> & A,
  const Eigen::DenseBase<DerivedB> & B,
  Eigen::PlainObjectBase<DerivedC> & C,
  Eigen::PlainObjectBase<DerivedIA> & IA,
  Eigen::PlainObjectBase<DerivedIB> & IB)
{
  DerivedC AB,BA;
  DerivedIA IAB,IBA;
  setdiff(A,B,AB,IAB);
  setdiff(B,A,BA,IBA);
  setunion(AB,BA,C,IA,IB);
  slice(IAB,DerivedIA(IA),IA);
  slice(IBA,DerivedIB(IB),IB);
}
